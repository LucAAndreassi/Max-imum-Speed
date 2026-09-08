# Project 1 - Design/Build/Fly Aircraft Design Optimization

## Problem Outline:
Design/Build/Fly is an annual collegiate aircraft design competition hosted by the AIAA, that several of our team members are participating in as part of an ASU Student Organization. The rules for this years' competition are linked [here](https://aiaa.org/wp-content/uploads/2026/09/DBF-2027-Rules-Draft.pdf)..

The objective is to design a remotely-piloted UAV to carry a 'sensor' payload for transportation and deployment. In our simplification of the ruleset, there are three scoring missions:

**Ground Mission**:  Drop the sensor within a protective shipping container from 5 feet off the ground and sustain no damage.
    GM score is the sensor weight in lbs, normalized to a maximum weight of 10 lbs

**Mission 2**: Carry the shipping container with sensor inside and fly 5 laps around the course within 5 minutes.
    M2 score is the sensor weight in lbs, divided by the time to complete the 5 laps in seconds. Normalized to 0.15 lb/s

**Mission 3**: Deploy the sensor(M3 sensor weight can be less than M2 and GM) in flight and fly as many laps as possible in 5 minutes.
    M3 score is the M3 sensor weight in lbs times the number of laps completed in 5 minutes. Normalized to 50 lb*laps

### Course:
The aircraft will take off from the 'Starting Line', and each lap is counted once the starting line is passed. Landing is not a part of the lap.

![Course Layout](image.png)
---

### Objective Function:
$$
 Score = \frac{SensorWeight_{1,2}}{10} + \frac{SensorWeight_{1,2}/Time}{0.15} + \frac{SensorWeight_{3}*Laps_3}{50}
$$
Success in the competition is purely determined by maximizing the total score, which is a sum of all three mission scores. This problem is non-convex, non-linear, and incorporates discrete variables such as # of laps.

---

### Constraints:

#### Rule-Imposed Limitations
- Wing span, $b \le 6 ft$
- Aircraft Weight, $W_T \le 55 lbs$
- Propulsion battery energy $\le 100Wh$

#### Operational Limitations
- Aircraft Weight $W_T \le 20 lbs$
- Aircraft Empty Weight $W_e \ge 10 lbs$

---

### Early Design Decisions

- Wing span will be 6 ft, to maximize span-wise efficiency
- Propulsion will be modelled after a known/tested motor and propeller combination.
- Aircraft configuration will be fixed-wing, with a single wing and inverted T-Tail
- Fuselage will have a circular cross-section
- Sensor will be a 6:1 length:diameter cylinder filled with leadshot with a nosecone for reduced drag
- Turns will be done at 2.5g
- Static margin at 10%
  
### Design Parameters to Optimize
- Wing Area, $S_w$
- Fuselage Length, $L_F$
- Sensor Weights, $SW_{1,2}$ and $SW_3$
- Cruise Velocities, $V_2$ and $V_3$

All other aircraft parameters, such as the sizing of the vertical/horizontal tail surfaces and the fuselage diameter, will be driven by these optimization variables.