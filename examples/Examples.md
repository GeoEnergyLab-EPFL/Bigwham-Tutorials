# Verification examples for the solution of Elastic crack problems with BigWham
### ROADMAP

- OK -> script exists and is in the repo already
- TBF -> some script exist somewhere (e.g. mathematica) but needs to be tidy up and finalized
- TBD -> yet to be done 

## 2D Tests

### Crack problems in an infinite medium 
+ Displacements and stresses around a unit dislocation segment [OK]
 + Griffith crack under remote tension
    - Displacement discontinuities from loading + convergence with 2DP0  [OK]  
    - Induced Stresses & displacements around an uniformly loaded crack  [OK]
    - Comparison 2DP0 and 2DP1 (convergence) [OK]  
    - unstructured mesh and robustness 2DP0 vs 2DP1 [OK]  
    - Simplified 3D P0 kernel test [OK] but not as jupytert
 + Dugdale-Barenblatt crack (remote tension)
    - displacement discontinuity from loading [TBF] - ratul
    - unstructured mesh and robustness 2DP0 vs 2DP1 [TBD]  
+ Arc crack under remote isotropic tension (mode I+II) - SIF  [TBD] - translate mma scripts
+ Star cracks under remote isotropic tension (mode I+II) - SIF [TBD] - translate mma scripts
+ Crack under a linear tensile stress gradient (Weertmaan's pulse) [TBD] - to be developed

### Elastic problems 
+ Circular cavity in an unbounded media under anisotropic far-field stress and internal pressure [TBF] - abhijeet 
+ Rectangle under tension  [TBF]
+ hollow cylinder with internal pressure [TBD] 
+ 'Complicated' geometry with mixed BC with the full traction BIE (comparison with FEM)  ? [TBD]

### Crack problems in bounded domain
+ Crack from a pressurized hole in an otherwise unbounded media - [TBF] abhijeet
+ Centered crack specimen (under remote uniaxial tension) [TBD]- tada part II 2.1 
+ double edge notch test specimen (Tada part II 2.6) [TBD]
+ Three point bending test (Tada) [TBD]


## 3D Tests

### Crack problems in an infinite medium
+ Penny-shaped crack under uniform remote loading (mode I or II/III)
    - Displacement discontinuities from loading (and vice-versa) 
        - mode I   with convergence     [TBD]
        - mode II/III  with convergence [TBD]
    - Comparisons between T6 and T0 (note T6 does not have the displacement) [TBD]
    - Stresses and displacement around mode I [OK] / II/III cracks [TBD]
+ Dugdale-Barenblatt penny-shaped crack under uniform remote mode I loading 
+ Penny-shaped crack under torsional loading 
        - comparisons T6 vs T0 [TBD]
+ Elliptical crack under remote tension [TBD]
+ Bowl-shaped crack under remote tension  [TBF]
+ 3D star cracks equivalent (long in one direction to compare with 2D) [TBD]

### Elastic problems 
+ A spherical cavity in an otherwise unbounded media ?   [TBD]
+ Square under tension (quarter model) with mixed Boundary Conditions [TBF]


### Crack problems in bounded domain
+ Penny-shaped crack under internal pressure in a finite cube  [TBD]
+ Penny-shaped crack parallel to a free-surface [TBD]
+ Penny-shaped crack intersecting a borehole [TBD]
+ Saw-cut cylindrical sample [TBD]


