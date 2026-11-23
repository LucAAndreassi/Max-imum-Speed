"""Own sparse Euler-Bernoulli bending FE and planar finite-rotation elastica.

SI throughout. Scalar bending subspace: [w, dw/dx] per node.
This is not a shell solver or a ply-level composite failure model.
"""
from dataclasses import dataclass
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu
from scipy.linalg import eigh
from vibration import first_frequency_partials, modal_summary, tip_mass_modes


def section(D, t):
    """Area and bending second moment; complex-safe for partial derivatives."""
    if np.real(D) <= 0 or np.real(t) <= 0 or np.real(2*t) >= np.real(D):
        raise ValueError("Require D > 2t > 0")
    d = D - 2*t
    return np.pi*(D*D-d*d)/4, np.pi*(D**4-d**4)/64


def hermite(xi, length):
    """Cubic interpolation on xi in [0,1], including nodal slopes."""
    return np.array([1-3*xi**2+2*xi**3,
                     length*(xi-2*xi**2+xi**3),
                     3*xi**2-2*xi**3, length*(-xi**2+xi**3)])


def element(length, EI, mass_per_length, q):
    """Consistent beam stiffness, translational inertia, and uniform load."""
    l = length
    k = EI/l**3*np.array([[12,6*l,-12,6*l], [6*l,4*l*l,-6*l,2*l*l],
                          [-12,-6*l,12,-6*l], [6*l,2*l*l,-6*l,4*l*l]])
    m = mass_per_length*l/420*np.array([
        [156,22*l,54,-13*l], [22*l,4*l*l,13*l,-3*l*l],
        [54,13*l,156,-22*l], [-13*l,-3*l*l,-22*l,4*l*l]])
    f = q*l/12*np.array([6,l,6,-l])
    return k, m, f


@dataclass
class Beam:
    length: float
    E: float
    mass_per_length: float
    tip_mass: float
    tip_force: float
    gravity: float = 9.80665
    elements: int = 20

    def matrices(self, D, t, free=None):
        _, I = section(D,t)
        n = self.elements
        nd = 2*(n+1)
        k, m, fe = element(self.length/n, self.E*I,
                           self.mass_per_length, -self.mass_per_length*self.gravity)
        rows, cols, ks, ms = [], [], [], []
        F = np.zeros(nd, dtype=np.result_type(k,fe))
        for e in range(n):
            ids = np.arange(2*e,2*e+4)
            rows.extend(np.repeat(ids,4)); cols.extend(np.tile(ids,4))
            ks.extend(k.ravel()); ms.extend(m.ravel())
            F[ids] += fe
        F[-2] += self.tip_force-self.tip_mass*self.gravity
        K = coo_matrix((ks,(rows,cols)),shape=(nd,nd)).tocsc()
        M = coo_matrix((ms,(rows,cols)),shape=(nd,nd)).tolil()
        M[-2,-2] += self.tip_mass
        if free is None:
            free = np.arange(2,nd)
        return K, M.tocsc(), F, np.asarray(free)

    def solve(self, D, t, free=None):
        K,M,F,free = self.matrices(D,t,free)
        lu = splu(K[free,:][:,free])
        u = np.zeros(len(F),dtype=F.dtype)
        u[free] = lu.solve(F[free])
        return u, K@u-F, (K,M,F,free,lu)

    def modes(self,D,t,count=4,free=None,method="naca_tn_3459"):
        """NACA TN 3459 equations (13)-(14), sampled at the FE nodes.

        Use method='fem' or fe_modes for the independent matrix eigensolve.
        The analytical method requires the standard clamped-free boundary.
        """
        if method == "fem":
            return self.fe_modes(D,t,count,free)
        if method != "naca_tn_3459":
            raise ValueError("Vibration method must be 'naca_tn_3459' or 'fem'")
        if free is not None:
            raise ValueError("NACA cantilever modes require the fixed-root boundary; use method='fem' for custom supports")
        _,I=section(D,t)
        return tip_mass_modes(self.length,self.E*I,self.mass_per_length,self.tip_mass,
                              np.linspace(0,self.length,self.elements+1),count)

    def vibration_summary(self,D,t,count=4):
        _,I=section(D,t)
        return modal_summary(self.length,self.E*I,self.mass_per_length,self.tip_mass,count)

    def fe_modes(self,D,t,count=4,free=None):
        """Independent consistent-mass FE eigenproblem, including tip mass."""
        K,M,_,free = self.matrices(D,t,free)
        vals, vecs = eigh(K[free,:][:,free].toarray(),M[free,:][:,free].toarray(),
                           subset_by_index=[0,min(count,len(free))-1])
        if np.min(vals) <= 0:
            raise ValueError("Nonpositive supported eigenvalue")
        modes = np.zeros((K.shape[0],len(vals)))
        modes[free] = vecs
        return np.sqrt(vals)/(2*np.pi), modes

    def end_moments(self,D,t,u):
        _,I = section(D,t)
        k,_,f = element(self.length/self.elements,self.E*I,
                         self.mass_per_length,-self.mass_per_length*self.gravity)
        return np.array([(k@u[2*e:2*e+4]-f)[[1,3]] for e in range(self.elements)])

    def tip_gradient(self,D,t):
        """One transposed solve, plus partial matrix derivatives; dm/dD=0.

        For a continuous density-based design use a mass law externally;
        stock product line mass is an independently specified measurement.
        """
        u,_,(K,M,F,free,lu) = self.solve(D,t)
        rhs = np.zeros(len(free)); rhs[-2] = 1
        adj = lu.solve(rhs,trans="T")
        grad = []
        for j in range(2):
            z = np.array([D,t],complex); z[j] += 1e-25j
            Kc,_,Fc,_ = self.matrices(*z)
            h_d = (Kc.imag@u-Fc.imag)/1e-25
            grad.append(-adj@h_d[free])
        return np.array(grad)

    def scalar_gradient(self,D,t,quantity):
        """Adjoint of a complex-safe quantity q(D,t,u), one state solve.

        Sign convention: K^T a = q_u, dq/dd = q_d - a^T(K_d u-F_d).
        Complex step evaluates partials, never replaces the state adjoint.
        """
        u,_,(_,_,_,free,lu)=self.solve(D,t)
        h=1e-25
        gu=np.zeros(len(free))
        for j,i in enumerate(free):
            uc=u.astype(complex);uc[i]+=1j*h
            gu[j]=np.imag(quantity(D,t,uc))/h
        adj=lu.solve(gu,trans="T")
        grad=[]
        for j in range(2):
            z=np.array([D,t],complex);z[j]+=1j*h
            K,_,F,_=self.matrices(*z)
            hd=(K.imag@u-F.imag)/h
            grad.append(np.imag(quantity(*z,u))/h-adj@hd[free])
        return np.array(grad)

    def frequency_gradient(self,D,t):
        """NACA implicit-root geometry gradient, holding stock line mass fixed."""
        _,I=section(D,t)
        _,partials=first_frequency_partials(self.length,self.E*I,self.mass_per_length,self.tip_mass)
        grad=[]
        for j in range(2):
            z=np.array([D,t],complex); z[j]+=1e-25j
            dEI=self.E*section(*z)[1].imag/1e-25
            grad.append(partials["EI"]*dEI)
        return np.array(grad)


class Elastica:
    """Inextensible planar beam, linear angle FE, four-point Gauss.

    Pi = int EI/2 theta'^2 ds - int W(s) sin(theta) ds - T theta(L)
    W(s) = tip force + integral_s^L q(v) dv (including tip weight).
    Optional follower force stays normal to the tip tangent. Its residual is
    assembled from virtual work, without an energy potential; tangent is
    nonsymmetric. tip_force and q remain fixed-direction dead forces.
    Positive vertical force rotates toward +z.
    """
    def __init__(self,length,EI,tip_force,q=0.,elements=40,tip_moment=0.,follower_force=0.):
        self.length,self.EI,self.P,self.q = length,EI,tip_force,q
        self.n,self.T=elements,tip_moment
        self.follower=follower_force
        x,w=np.polynomial.legendre.leggauss(4)
        self.x=(x+1)/2; self.w=w/2

    def residual(self,theta,EI=None):
        if EI is None: EI=self.EI
        n=self.n; l=self.length/n
        R=np.zeros(n+1,dtype=np.result_type(theta,EI))
        H=np.zeros((n+1,n+1),dtype=R.dtype)
        energy=0.*EI
        B=np.array([-1.,1.])/l
        for e in range(n):
            ids=np.array([e,e+1]); a=theta[ids]
            re=EI*l*(B@a)*B
            he=np.asarray(EI*l*np.outer(B,B),dtype=R.dtype)
            energy+=EI*l*(B@a)**2/2
            for x,w in zip(self.x,self.w):
                N=np.array([1-x,x]); th=N@a
                W=self.P+self.q*(self.length-(e+x)*l)
                delta=theta[-1]-th
                re-=l*w*(W*np.cos(th)+self.follower*np.cos(delta))*N
                he+=l*w*W*np.sin(th)*np.outer(N,N)
                he-=l*w*self.follower*np.sin(delta)*np.outer(N,N)
                H[ids,-1]+=l*w*self.follower*np.sin(delta)*N
                energy-=l*w*W*np.sin(th)
            R[ids]+=re; H[np.ix_(ids,ids)]+=he
        R[-1]-=self.T; energy-=self.T*theta[-1]
        return R,H,energy

    def solve(self,steps=10,tolerance=1e-11):
        theta=np.zeros(self.n+1)
        original=(self.P,self.q,self.T,self.follower)
        records=[]
        try:
            for load in np.linspace(1/steps,1,steps):
                self.P,self.q,self.T,self.follower=[v*load for v in original]
                residuals=[]
                for iteration in range(40):
                    R,H,energy=self.residual(theta)
                    norm=float(np.linalg.norm(R[1:],np.inf));residuals.append(norm)
                    scale=max(1.,abs(self.P)*self.length,abs(self.follower)*self.length,abs(self.q)*self.length**2,abs(self.T))
                    if norm<=tolerance*scale: break
                    step=np.linalg.solve(H[1:,1:],-R[1:])
                    descent=R[1:]@step
                    a=1.
                    for _ in range(30):
                        trial=theta.copy();trial[1:]+=a*step
                        Rt,_,et=self.residual(trial)
                        # Residual decrease also accepts final roundoff-limited steps.
                        if (not self.follower and et<=energy+1e-4*a*descent) or np.linalg.norm(Rt[1:],np.inf)<(1-1e-4*a)*norm:
                            theta=trial;break
                        a*=0.5
                    else: raise RuntimeError("Elastica line search failed")
                else: raise RuntimeError("Elastica Newton did not converge")
                records.append({"load_factor":float(load),"residuals":residuals})
        finally:
            self.P,self.q,self.T,self.follower=original
        if not self.follower and np.min(np.linalg.eigvalsh(self.residual(theta)[1][1:,1:]))<=0:
            raise RuntimeError("Elastica tangent is not positive definite")
        return theta,records

    def coordinates(self,theta):
        out=[np.zeros(2,dtype=np.result_type(theta))]
        l=self.length/self.n
        for e in range(self.n):
            v=np.zeros(2,dtype=np.result_type(theta))
            for x,w in zip(self.x,self.w):
                th=(1-x)*theta[e]+x*theta[e+1]
                v+=l*w*np.array([np.cos(th),np.sin(th)])
            out.append(out[-1]+v)
        return np.array(out)

    def tip_EI_gradient(self,theta):
        """Implicit adjoint d z_tip / d EI, with converged angle held in partials."""
        h=1e-25
        gu=np.zeros(self.n)
        for i in range(1,self.n+1):
            tc=theta.astype(complex);tc[i]+=1j*h
            gu[i-1]=self.coordinates(tc)[-1,1].imag/h
        _,H,_=self.residual(theta)
        adj=np.linalg.solve(H[1:,1:].T,gu)
        Rd=self.residual(theta,self.EI+1j*h)[0].imag/h
        return -adj@Rd[1:]

    def coordinates_at(self,theta,samples):
        """Integrate the element's linear angle field to requested arc lengths."""
        nodes=self.coordinates(theta)
        l=self.length/self.n
        out=[]
        for s in samples:
            if s<0 or s>self.length:raise ValueError("Sample outside beam")
            e=min(int(s/l),self.n-1);fraction=(s-e*l)/l
            xy=nodes[e].copy()
            for x,w in zip(self.x,self.w):
                xi=x*fraction
                th=(1-xi)*theta[e]+xi*theta[e+1]
                xy+=l*fraction*w*np.array([np.cos(th),np.sin(th)])
            out.append(xy)
        return np.array(out)
