# Research Review: Dr. Kiran Ramesh

**Associate Professor of Aerospace Engineering**
School for Engineering of Matter, Transport and Energy (SEMTE), Arizona State University, Tempe

> Prepared for Project 3 team reference. Compiled 18 October 2026.

---

## Contents

- [Plain-language summary](#plain-language-summary)
- [Paper 1 — LESP and LDVM (2014)](#paper-1--lesp-and-ldvm-2014)
- [Paper 2 — Morphing airfoils (2024)](#paper-2--morphing-airfoils-2024)
- [Paper 3 — Wake vorticity and the Kutta condition (2025)](#paper-3--wake-vorticity-and-the-kutta-condition-2025)
- [The through-line](#the-through-line)
- [What we could not find](#what-we-could-not-find)
- [Field context](#field-context)
- [References](#references)

---

## Plain-language summary

### The problem he works on

When a wing moves gently through steady air, predicting lift is a solved problem — you can do it with equations from the 1930s. But when a wing pitches sharply, flaps, changes shape, or flies into a gust, the air separates off the front edge and rolls up into a vortex. That vortex sits on top of the wing and changes the forces dramatically: sometimes doubling the lift for a moment, then dumping it when the vortex washes away. This is what insects, birds and bats exploit to fly, and it is what makes small drones, helicopter rotors and wind-turbine blades hard to design.

The honest way to predict this is full CFD, which takes hours or days per run — far too slow for a design loop where you want to try thousands of shapes or motions. The field needs something in between: fast enough to run in a loop, accurate enough to trust.

### What he built

Ramesh's career contribution is exactly that middle tool. The headline idea is a single number, the **Leading-Edge Suction Parameter (LESP)**, measuring how hard the air is being pulled around the front edge of the wing. The insight: when LESP crosses a critical threshold a vortex starts forming, and below it the flow stays attached. That one number converts a messy, hard-to-predict event into a simple on/off criterion checkable at every time step.

He wrapped that criterion into a method (**LDVM**) representing the wing and its wake as a cloud of discrete vortex points, releasing a new one whenever LESP goes critical. It runs in tens of seconds instead of hours and matches CFD and wind-tunnel experiments well.

### Where he is now

| Direction | Year | What it is |
|---|---|---|
| **Morphing wings** | 2024 | Change the wing's shape mid-flight to control the vortex on purpose — make it form earlier, later, stronger, weaker. A design and control question, and our most relevant thread. |
| **Getting the math right** | 2025 | A technical correction to how the shed wake is represented each time step, making the model both more rigorous and cheaper. |

### Why this matters for our project

> [!IMPORTANT]
> **He already has a fast surrogate.** We cannot walk in and offer him speed. Anything we propose has to address either accuracy where his model is currently wrong, or a search problem his fast model still cannot make tractable.

---

## Paper 1 — LESP and LDVM (2014)

> Ramesh, K., Gopalarathnam, A., Granlund, K., Ol, M. V. & Edwards, J. R. (2014). "Discrete-vortex method with novel shedding criterion for unsteady aerofoil flows with intermittent leading-edge vortex shedding." *Journal of Fluid Mechanics*, **751**, 500–538.

### What problem it solves

Existing low-order methods either assumed vortices shed continuously (valid only for sharp leading edges) or used ad-hoc criteria like a fixed critical angle of attack (valid only for a narrow set of motions). Classical unsteady airfoil theories were limited to small amplitudes and attached leading-edge flow. Neither covered the realistic case: a rounded leading edge with vortices forming and stopping intermittently.

### Key technical points

- Introduces **LESP** as a criterion for intermittent shedding from rounded leading edges, governed by a maximum allowable leading-edge suction.
- Shows that within unsteady thin-airfoil theory, LESP is the $A_0$ term of the Fourier series representing the chordwise variation of bound vorticity. *This is the elegant part — the criterion isn't bolted on, it falls out of the existing theory.*
- When instantaneous LESP exceeds the critical value, a discrete vortex is released at the leading edge, with strength set so LESP is held **at** the critical value during shedding.
- A single empirical parameter (critical LESP) determines onset, growth and termination of the leading-edge vortex.
- Validated against CFD and experiment across motions from the unsteady airfoil literature.

### The open soft spot

Critical LESP is empirical. A companion 2012 AIAA paper (Ramesh, Ke, Gopalarathnam & Edwards) showed it depends on **airfoil shape and Reynolds number** — increasing leading-edge radius raises the critical value, because a more rounded edge supports more suction before separating. They tabulated values across a parameter range.

> [!NOTE]
> **Confirm in the interview** whether this is still a table lookup from CFD.

### Related follow-on

SureshBabu, Ramesh & Gopalarathnam (2019), *AIAA Journal* **57**(4), 1409–1422, "Model reduction in discrete-vortex methods for unsteady airfoil flows." Merges pairs of discrete vortices under a condition preserving the velocity at the leading edge, conserving bound circulation and leading-edge suction. Same accuracy, significantly lower runtime.

**Why it matters to us:** he has *already* optimized his own method's cost.

---

## Paper 2 — Morphing airfoils (2024)

> Martínez-Carmena, A. & Ramesh, K. (2024). "Inviscid modeling of unsteady morphing airfoils using a discrete-vortex method." *Theoretical and Computational Fluid Dynamics*, **38**(6), 845–862. [doi:10.1007/s00162-023-00678-7](https://doi.org/10.1007/s00162-023-00678-7) — open access, EPSRC-funded.

### What problem it solves

Thin-airfoil theory assumes small flow disturbances, which breaks as soon as the wing itself deforms substantially. Morphing wings — variable camber, trailing-edge flaps, anything bio-inspired that bends — fall outside that assumption. This paper extends the low-order model to large-amplitude camber change.

### Key technical points

- Combines potential-flow theory adapted for unsteady airfoils with discrete-vortex numerical methods, for rapid predictions of flow behavior and force evolution.
- **Core trick:** a *time-varying chord line* over which the boundary condition is satisfied, rather than a fixed one. This is what escapes thin-airfoil theory's small-disturbance restriction and allows large camber-line deformation.
- The difference from the fixed-chord formulation is most pronounced for large-amplitude deformations.
- Author keywords: morphing airfoil, variable camber, trailing-edge flap, vortex methods, thin-airfoil theory, bio-inspired.

### The design angle — why this is our best thread

The associated PhD thesis (Martínez Carmena, *Vortex flows around morphing foils: reshape yourself to effect a change*, University of Glasgow 2023, supervised by Ramesh) states the goal directly: demonstrate that camber morphing can affect leading-edge vortices, with the amount of vorticity produced at the leading edge modulated by suitably tailoring the shape of the camber line. The thesis also notes that the strength of nascent particles at the leading edge is computed from the velocity at the edge of the shear layer, formulated in terms of LESP.

In other words: **shape parameters in, vortex characteristics out.** That is a parametric design map, which is a many-query problem, which is what our class is about.

> [!TIP]
> Lead the interview here. Ask how he searches that design space today and how many model runs one search takes.

---

## Paper 3 — Wake vorticity and the Kutta condition (2025)

> Ramesh, K. (2025). "Correct Wake Vorticity Shedding in Time-Stepping Unsteady Thin-Airfoil Theory." *AIAA Journal*. [doi:10.2514/1.J065056](https://doi.org/10.2514/1.J065056) — sole-authored, listed under Arizona State University.

> [!WARNING]
> **Sourcing caveat.** This paper is paywalled; we have only the title, affiliation and subject tags. The content below comes from its free precursor preprint, [arXiv:2205.08647](https://arxiv.org/abs/2205.08647), *"On Satisfying the Kutta Condition in Unsteady Thin Aerofoil Theory"* (May 2022), covering the same contribution. Treat as closely related, not identical — and say so to him rather than implying we read it.

### What problem it solves

In unsteady thin-airfoil theory you must apply a Kutta condition at the trailing edge to pin down the circulation around the airfoil. How that is done in a time-stepping scheme turns out to matter, and the conventional treatment is not strictly correct.

### Key technical points

- Conventional approach represents the wake shed at each discrete time step as a point vortex or regularized vortex blob.
- His approach instead represents it as a **continuous distribution of vorticity derived from the exact Wagner solution.**
- Demonstrated on three canonical cases: step change in angle of attack (Wagner problem), harmonic heaving (Theodorsen problem), and a pitch-ramp-hold maneuver.
- Verified against exact analytical solutions *and* CFD of the incompressible Euler equations.
- Result: Kutta condition satisfied at all reduced frequencies, with finite velocities and pressure difference going to zero at the trailing edge.
- Claims improvement on three axes at once — theoretical rigor, computational cost, and numerical accuracy.

### What it tells us about him

He is actively working on the **correctness** of low-order methods, not just their speed. A pitch framed as "we made it faster" will land badly. "We found where it is wrong" is the better frame.

### Useful side benefit

The Wagner and Theodorsen problems have **exact analytical solutions**. If we need a bulletproof reference solution for rubric item P1, those are free.

---

## The through-line

Across all three papers the pattern is identical:

```
START      classical potential-flow theory (Wagner, Theodorsen, thin-airfoil)
  ↓
ADD        discrete-vortex elements to capture shed vorticity
  ↓
FIX        the one place the classical theory breaks
  ↓
VALIDATE   against CFD and experiment
  ↓
GOAL       a tool fast and rigorous enough for design and control loops
```

What gets fixed, paper by paper:

| Year | What classical theory couldn't do | His fix |
|---|---|---|
| 2014 | Handle intermittent leading-edge shedding | LESP criterion |
| 2024 | Handle large shape change | Time-varying chord line |
| 2025 | Represent the wake correctly at the trailing edge | Wagner-based continuous vorticity distribution |

He is not a CFD person and not an ML person. He is a **theory-plus-low-order-modeling** person whose stated ambition at ASU includes artificial intelligence, but who has — as far as we can find — no published ML work.

---

## What we could not find

Be honest about these gaps in the report.

- **No published ML or neural-network work** by him that we could locate. His SciML signal is best described as *Adjacent* (stated AI ambition, no output yet), not *Active*.
- **No ASU lab website or current student roster.** He is new here (2024–25), so this may simply not exist yet.
- **No full text** of the 2025 *AIAA Journal* paper.
- **No public statement of his computational bottleneck.** This is the central thing the interview needs to establish.

---

## Field context

"Extreme aerodynamics: A data-driven perspective," *Physical Review Fluids* **11**, 014702, published 15 January 2026, cites his 2014 LDVM paper among its references. The review frames the problem as small air vehicles operating in urban canyons, mountainous terrain and ship wakes, where gust scales can exceed the vehicle's own aerodynamic scales.

**Takeaway:** data-driven methods for exactly his class of problem are being actively surveyed as of nine months ago. He will be aware of this. Don't explain the landscape to him — ask where he thinks it's overhyped.

---

## References

| Source | Link |
|---|---|
| ASU profile | https://search.asu.edu/profile/1348363 |
| ASU welcome piece | https://fullcircle.asu.edu/welcome/kiran-ramesh |
| 2014 — LESP / LDVM | https://eprints.gla.ac.uk/99370 |
| 2019 — Model reduction | https://eprints.gla.ac.uk/183312 |
| 2024 — Morphing airfoils | https://doi.org/10.1007/s00162-023-00678-7 |
| 2025 — AIAA Journal | https://doi.org/10.2514/1.J065056 |
| 2022 — Free preprint | https://arxiv.org/abs/2205.08647 |
| 2023 — Student thesis | https://theses.gla.ac.uk/83534 |
| 2026 — PRF review | https://journals.aps.org/prfluids/abstract/10.1103/8pr4-txmt |
