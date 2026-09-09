#This Class
import numpy as np
import math

class PlaneModel:
    def __init__(self, wing_area, fuselage_length,SW2):
        self.wing_area = wing_area  # in square feet
        self.fuselage_length = fuselage_length  # in feet
        self.SW2 = SW2  # in lbs

        # Sizing based on input parameters, determined by helper functions below
        self.sensor_size()
        self.geo_sizing()


    def sensor_size(self):
        """
        Calculate the sensor size based on the maximum weight and given aspect ratio.
        :return: Sensor diameter in square feet, and the total drag area
        """
        aspect_ratio = 6  # Length to width ratio of the sensor
        leadshot_density = 430 # lbs/ft^3
        
        sensor_volume = self.SW2 / leadshot_density
        sensor_length = (sensor_volume*(aspect_ratio*2)**2/math.pi) ** (1 / 3)  # Use given aspect ratio to find length
        
        Cd0nose = 0.2 # Zero-lift drag coefficient for the nose
        Cfbody = 0.0042 # Skin friction coefficient for the cylinder body

        # Outputs
        self.sensor_diameter = sensor_length / aspect_ratio  # Calculate diameter based on length and aspect ratio
        self.sensor_drag_area = Cd0nose * self.sensor_diameter**2 * math.pi / 4  + Cfbody*math.pi*self.sensor_diameter*sensor_length # Calculate the total drag area
    
    def geo_sizing(self):
        """
        Calculate the size of aero surfaces, and fuselage. Used to determine weight and drag of the aircraft.
    
        """
        fuselage_diameter = self.sensor_diameter*1.5 # Buffer size so the sensor fits inside the fuselage
        self.aspect_ratio = 6**2 / self.wing_area
        self.e = 0.75 # Assume rectangle, can mess around with increasing span efficiency in detailed design
        self.Kwing = 1/math.pi/self.e/self.aspect_ratio
        self.chord = 6 / self.aspect_ratio

        # Mass props
        fuselage_weight = 3.5*math.pi*(fuselage_diameter/2)**2*self.fuselage_length # Based on last year's competition sizing
        electronics_weight = 2.5 # Based on last year's competition sizing
        wing_weight = 2/7*self.wing_area # Based on last year's competition sizing
        empennage_weight = 1 # Based on last year's competition sizing
        
        self.total_weight = fuselage_weight + electronics_weight + wing_weight + empennage_weight
        self.cgx = (electronics_weight*0.5 + empennage_weight*self.fuselage_length + fuselage_weight*self.fuselage_length/2 + wing_weight*self.chord/10) / (self.total_weight - wing_weight) # Assume wing is at 1/10 chord, electronics at 0.5 ft, empennage at end of fuselage, and fuselage weight at 1/2 fuselage length
        print(self.total_weight)
        # Tail Sizing
        lt = self.fuselage_length - (self.cgx + self.chord/10) # Distance from aerodynamic center to tail
        self.Sv = 0.07*self.wing_area*6/lt 
        self.Sh = 0.7*self.wing_area*self.chord/lt

        # Drag Area for everything but the wing
        cd0fuselage = 0.3 # Nose and propeller frontal drag coefficient
        cfskin = 0.0042 # Skin friction coefficient for all external surfaces
        self.drag_area = cd0fuselage * fuselage_diameter**2 * math.pi / 4 + cfskin * (math.pi*fuselage_diameter*self.fuselage_length + 2*(self.Sh + self.Sv)) # Drag area for fuselage and empennage



    def calculate_power(self,velocity, lift, sensor_deployed):
        # First need to calculate drag, which we need the wing lift and ultimately velocity for dynamic pressure q

        rho = 0.0023769 # Slug/ft^3, density of air at sea level
        q = 0.5 * rho * velocity**2 # Dynamic pressure

        cl = lift / (q * self.wing_area) # Lift coefficient
        clmax = 1.2
        if cl > clmax:
            cl = clmax
            return False
        cd_induced = self.Kwing * cl**2 # Induced drag coefficient
        cd_wing = 0.02 + cd_induced # Total drag coefficient for the wing
        wing_drag = cd_wing * q * self.wing_area # Drag force on the wing
        body_drag = self.drag_area * q # Drag force on the fuselage and empennage
        if sensor_deployed == True:
            sensor_drag = self.sensor_drag_area * q # Drag force on the sensor
        else:
            sensor_drag = 0

        total_drag = wing_drag + body_drag + sensor_drag # Total drag force on the aircraft

        # Now that we have drag, we can calculate the power required to generate the thrust at that given velocity
        th = np.linspace(0.3, 1, 100)
        T = 21.175 
        V = velocity 
        Do = total_drag/velocity**2
        Vmax = 170 # pitch speed in ft/s
        v = (-T * (1.4 * th - 0.4) / (Do * Vmax *th) + np.sqrt((T * (1.4 * th - 0.4) / (Do * Vmax * th))**2 + 4 * T * (1.4 * th - 0.4) / Do)) / 2
        throttle = th[np.abs(v-V) == np.min(np.abs(v-V))]
        Thrust = T*(1.4*throttle-0.4)
        Amps = 85*(1.8541 * throttle**2 - 1.0736 * throttle + 0.2136)
        
        return Amps



class Course:

    def __init__(self, wing_area, fuselage_length, SW2, SW3, V2, V3):
        self.wing_area = wing_area
        self.fuselage_length = fuselage_length
        self.SW2 = SW2
        self.SW3 = SW2 - abs(SW3)
        self.V2 = V2
        self.V3 = V3

    def calculate_forces(self):
        aircraft = PlaneModel(self.wing_area, self.fuselage_length, self.SW2)
        unloaded_weight = aircraft.total_weight
        M2_weight = unloaded_weight + self.SW2*1.1 # Accounting for shipping container weight as a function of sensor weight
        M3_weight = unloaded_weight + self.SW3 

        LwingM2 = M2_weight - (aircraft.cgx + aircraft.chord/10) * M2_weight / (aircraft.cgx + aircraft.chord/10 - self.fuselage_length) # Lift on the wing at M2
        LwingM3 = M3_weight - (aircraft.cgx + aircraft.chord/10) * M3_weight / (aircraft.cgx + aircraft.chord/10 - self.fuselage_length) # Lift on the wing at M3
        AmpsM2 = aircraft.calculate_power(self.V2, LwingM2, sensor_deployed=False) # Amps at M2
        AmpsM3 = aircraft.calculate_power(self.V3, LwingM3, sensor_deployed=True) # Amps at M3

        # For turns, assume 2.5g turn
        self.nM2 = 2.5
        self.nM3 = 2.5
        LwingM2_turn = self.nM2 * LwingM2
        LwingM3_turn = self.nM3 * LwingM3
        AmpsM2_turn = aircraft.calculate_power(self.V2, LwingM2_turn, sensor_deployed=False) # Amps at M2 turn
        if AmpsM2_turn == False:
            result = False
            while result == False:
                self.nM2 = self.nM2 - 0.1
                LwingM2_turn = self.nM2 * LwingM2
                result = aircraft.calculate_power(self.V2, LwingM2_turn, sensor_deployed=False)
            AmpsM2_turn = result
        AmpsM3_turn = aircraft.calculate_power(self.V3, LwingM3_turn, sensor_deployed=True) # Amps at M3 turn
        if AmpsM3_turn == False:
            result = False
            while result == False:
                self.nM3 = self.nM3 - 0.1
                LwingM3_turn = self.nM3 * LwingM3
                result = aircraft.calculate_power(self.V3, LwingM3_turn, sensor_deployed=True)
            AmpsM3_turn = result

        return AmpsM2, AmpsM3, AmpsM2_turn, AmpsM3_turn

    def calculate_lap_consumption(self, AmpsM2, AmpsM3, AmpsM2_turn, AmpsM3_turn):
        # Calculate the number of laps that can be completed based on the battery capacity and the current draw at each segment of the course
        
        straight_timeM2 = 2000 / self.V2  # Time to complete the straight segment at M2 in seconds
        straight_timeM3 = 2000 / self.V3  # Time to complete the straight segment at M3 in seconds
        turn_radiusM2 = self.V2**2 / (32.2*(self.nM2**2-1))
        turn_radiusM3 = self.V3**2 / (32.2*(self.nM3**2-1))
        
        turn_timeM2 = 4*math.pi*turn_radiusM2 / self.V2  # Time to complete all turns in M2 in seconds
        turn_timeM3 = 4*math.pi*turn_radiusM3 / self.V3  # Time to complete all turns in M3 in seconds

        total_capacity_drawM2 = (AmpsM2*(straight_timeM2) + AmpsM2_turn*(turn_timeM2))  # Average current draw over the course
        total_capacity_drawM3 = (AmpsM3*(straight_timeM3) + AmpsM3_turn*(turn_timeM3))  # Average current draw over the course
        
        lap_timeM2 = straight_timeM2 + turn_timeM2  # Average speed over the course
        lap_timeM3 = straight_timeM3 + turn_timeM3  # Average speed over the course
        print(lap_timeM2,lap_timeM3)
        return total_capacity_drawM2, total_capacity_drawM3, lap_timeM2, lap_timeM3

    def mission_scores(self):
        battery_capacity = 3300*0.75*3600/1000  # in mAh, and accounting for a 5% buffer for takeoff and landing

        # M2, where 5 laps have to be completed in 5 minutes
        
        raw_forces = self.calculate_forces()
        AmpsM2, AmpsM3, AmpsM2_turn, AmpsM3_turn = raw_forces
        total_capacity_drawM2, total_capacity_drawM3, lap_timeM2, lap_timeM3 = self.calculate_lap_consumption(AmpsM2, AmpsM3, AmpsM2_turn, AmpsM3_turn)
        M2_time = 5*lap_timeM2
        M2_capacity = total_capacity_drawM2*5
        '''
        if M2_time > 300:
            return np.nan, np.nan, np.nan
        if M2_capacity > battery_capacity:
            return np.nan, np.nan, np.nan
        '''
        M2_score = self.SW2/M2_time

        # M3, where 5 minutes are given to complete as many laps as possible

        num_lapsM3 = math.floor(300 / lap_timeM3)
        feasible_lapsM3 = np.floor(battery_capacity / total_capacity_drawM3)
        if num_lapsM3 > feasible_lapsM3:
            num_lapsM3 = feasible_lapsM3

        M3_score = self.SW3*num_lapsM3

        GM_score = self.SW2 

        return M2_score, M3_score, GM_score, num_lapsM3




def objective_function(params):
    wing_area, fuselage_length, SW2, SW3, V2, V3 = params
    
    course = Course(wing_area, fuselage_length, SW2, SW3, V2, V3)
    scores = course.mission_scores()
    M2_score, M3_score, GM_score, num_lapsM3 = course.mission_scores()
    
    score = - (M2_score/0.15 + M3_score/50 + GM_score/10)
    print(f"Testing -> Area: {wing_area:.2f}, Length: {fuselage_length:.2f}, SW2: {SW2:.2f}, SW3: {SW3:.2f}, V2: {V2:.1f}, V3: {V3:.1f}, M3 Laps: {num_lapsM3}, Score = {score.item():.2f}")
    return score
