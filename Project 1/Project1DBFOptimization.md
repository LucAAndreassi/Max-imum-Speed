# Project 1 - Design/Build/Fly Aircraft Design Optimization

## Problem Outline:
Design/Build/Fly (DBF) is an annual collegiate aircraft design competition hosted by the AIAA (American Institute of Aeronautics and Astronautics), that several of our team members are participating in as part of an ASU Student Organization. The rules for this years' competition are linked [here](https://aiaa.org/wp-content/uploads/2026/09/DBF-2027-Rules-Draft.pdf). The purpose of this event is to put teams through a full design cycle from initial requirements to a fully-built and mission capable system in the span of an academic year. A significant part of the design cycle is optimization to maximize scoring potential, and thus placement in the competition. Aside from the competition itself, there are many implications for the use of design optimizaiton within the UAV (Unmanned Aerial Vehicle) domain that share similar processes and applications with DBF. 

The objective for this year's competition is to design a remotely-piloted UAV to carry a 'sensor' payload for transportation and deployment. In our simplification of the ruleset, there are three scoring missions:

**Ground Mission**:  Drop the sensor within a protective shipping container from 5 feet off the ground and sustain no damage. This is done in a stationary position on the ground. 
    GM score is the sensor weight in lbs, normalized to a maximum weight of 10 lbs

**Mission 2**: Carry the shipping container with sensor inside and fly 5 laps around the course within 5 minutes.
    M2 score is the sensor weight in lbs, divided by the time to complete the 5 laps in seconds. Normalized to 0.15 lb/s

**Mission 3**: Deploy the sensor (M3 sensor weight can be less than M2 and GM) in flight and fly as many laps as possible in 5 minutes.
    M3 score is the M3 sensor weight in lbs times the number of laps completed in 5 minutes. Normalized to 50 lb*laps

### Course:
The aircraft will take off from the 'Starting Line', and each lap is counted once the starting line is passed. Landing is not a part of the lap. The course is an oval shape with a 360 degree turn on the straight opposite of the starting position. 

![Course Layout](image.png)
---

### Decision Variables:
- Wing Area, $S_w$
- Fuselage Length, $L_F$
- Sensor Weights, $SW_{1,2}$ and $SW_3$
- Cruise Velocities, $V_2$ and $V_3$

All other aircraft parameters, such as the sizing of the vertical/horizontal tail surfaces and the fuselage diameter, will be driven by these optimization variables.

### Objective Function:
$$
 Score = \frac{SensorWeight_{1,2}}{10} + \frac{SensorWeight_{1,2}/Time}{0.15} + \frac{SensorWeight_{3}*Laps_3}{50}
$$
Success in the competition is purely determined by maximizing the total score, which is a sum of all three mission scores. 

### Constraints:

#### Rule-Imposed Limitations
- Wing span, $b \le 6 ft$
- Aircraft Weight, $W_T \le 55 lbs$
- Propulsion battery energy $\le 100Wh$

#### Operational Limitations and Early Design Decisions 
These are design decisions based on engineering intuition and historically successful DBF teams. This is intended to reduce the computational load of the model, and enforce several non-optimization related parameters such as manufacturability and available hardware.   
- Aircraft Weight, $W_T \le 20 lbs$
- Aircraft Empty Weight, $W_e \ge 10 lbs$
- Wing span, $b = 6 ft$

#### Lap Simulator Simplifications
These are simplifications to estimate preliminary design and increase fidelity in certain parts of the simulation. 
- Propulsion will be modelled after a known/tested motor and propeller combination.
- Aircraft configuration will be fixed-wing, with a single wing and inverted T-Tail
- Fuselage will have a circular cross-section
- Sensor will be a 6:1 length:diameter cylinder filled with leadshot with a nosecone for reduced drag
- Turns will be done at 2.5g
- Static margin at 10%
  

### Problem Classification
This problem is non-convex, non-linear, and mixed integer. Within the scope of our design space, this problem is non-convex due to the uncertainty that only one solution exists based on the limited design space. If the design space is expanded to include all the variables, then it is possible to show the design is the global minimum (convex), but since it is not in the scope of this problem, this remains classified as a non-convex problem. Additionally, this problem is non-linear due to the nature of the lap simulator used to model the vehicle's performance as a function of its characteristic parameters. Various parts of the simulation such as lift and drag are non-linear by nature and result in this problem being non-linear. Finally, this problem is MINLP since various input parameters can only be expressed in integer form, such as the number of laps, and other values are allowed to be continuous, such as wing area, cruise speed, coefficients of lift and drag.  

### Solution Methodology
The design optimization is based off the results of a lap simulator intended to model the performance of the vehicle under certain parameters. The optimization itself varies the decision variables and imposes contraints, finally using the objective function to evaluate a configuration.  

### Results and Interpretation

### Code and Reproducibility