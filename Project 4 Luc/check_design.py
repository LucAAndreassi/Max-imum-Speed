"""Editable material screen. Never treats missing properties as passing."""
import argparse
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parent


def evaluate_strength(stress_Pa,shear_Pa,material,required_sf=2.):
    if required_sf<=0 or stress_Pa<0 or shear_Pa<0:
        raise ValueError("Require positive safety factor and nonnegative stress magnitudes")
    requirements=[("tension","design_tensile_strength_Pa",stress_Pa),
                  ("compression","design_compressive_strength_Pa",stress_Pa),
                  ("transverse_shear","design_shear_strength_Pa",shear_Pa)]
    checks={}
    for name,key,demand in requirements:
        strength=material.get(key)
        if strength is None:
            checks[name]={"reserve_factor":None,"status":"UNRESOLVED","reason":f"Set material.{key}"}
        elif strength<=0:
            raise ValueError(f"{key} must be positive or null")
        else:
            reserve=strength/demand if demand else None
            checks[name]={"reserve_factor":reserve,"status":"PASS" if not demand or reserve>=required_sf else "FAIL"}
    valid=[x["reserve_factor"] for x in checks.values() if x["reserve_factor"] is not None]
    checks["nominal_maximum_stress_screen"]={
        "status":"FAIL" if any(x["status"]=="FAIL" for x in checks.values()) else
                 "UNRESOLVED" if any(x["status"]=="UNRESOLVED" for x in checks.values()) else "PASS",
        "minimum_known_reserve_factor":min(valid) if valid else None}
    checks["overall_structure"]={"status":"UNRESOLVED","reason":"Local shell failure, joints, laminate interaction and dynamic stability not evaluated"}
    return checks


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs",type=Path,default=ROOT/"inputs.json")
    parser.add_argument("--results",type=Path,default=ROOT/"assets"/"results.json")
    args=parser.parse_args()
    inputs=json.loads(args.inputs.read_text())
    results=json.loads(args.results.read_text())
    if results["inputs"]!=inputs:
        raise SystemExit("Inputs changed: regenerate results with python run.py before checking")
    output={c["name"]:evaluate_strength(c["peak_bending_stress_Pa"],c["homogeneous_shear_screen_Pa"],
                                       inputs["material"],inputs["required_safety_factor"])
            for c in results["candidates"]}
    print(json.dumps(output,indent=2,allow_nan=False))
    return 2 if any(c["overall_structure"]["status"]=="UNRESOLVED" for c in output.values()) else 0


if __name__=="__main__":raise SystemExit(main())
