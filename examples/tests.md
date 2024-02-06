# Verification tests for Elastic crack problems



## 2D Tests

### Crack problems in an infinite medium 
+ Displacements and stresses around a unit dislocation segment [OK]
 + Griffith crack under remote tension
    - Displacement discontinuities from loading   [OK] - see ratul
    - loading from displacement discontinuities   [OK] - see ratul
    - Induced Stresses around an uniformly loaded crack  [OK]
    - Comparison 2DP0 and 2DP1 (convergence) [OK] - see ratul
 + Dugdale-Barenblatt crack (remote tension)
    - displacement discontinuity from loading [ ] - see ratul 
+ Arc crack under remote isotropic tension (mode I+II) - SIF   - translate mma scripts
+ Star cracks under remote isotropic tension (mode I+II) - SIF - translate mma scripts
+ Crack under a linear tensile stress gradient (Weertmaan's pulse) - to be developed

### Elastic problems 
+ Circular cavity in an unbounded media under anisotropic far-field stress and internal pressure - abhijeet
+ Square under tension (quarter model vs full-model) 
+ hollow cylinder with internal pressure 
+ 'Complicated' geometry with mixed BC (comparison with FEM)  ?

### Crack problems in bounded domain
+ Crack from a pressurized hole in an otherwise unbounded media - abhijeet
+ Centered crack specimen (under remote uniaxial tension) - tada part II 2.1 
+ double edge notch test specimen (Tada part II 2.6)
+ Three point bending test (Tada)


## 3D Tests

### Crack problems in an infinite medium
+ Penny-shaped crack under remote loading (mode I or II/III)
    - Displacement discontinuities from loading (and vice-versa) (mode I or II/III)
    - Stresses and displacement around mode I [ok] / II cracks [tbd]
    - Comparisons between T6 and T0 (note T6 does not have the displacement) 
+ Dugdale-Barenblatt penny-shaped crack under remote mode I loading 
+ Elliptical crack under remote tension 
+ Bowl-shaped crack 

### Elastic problems 
+ A spherical cavity in an otherwise unbounded media ?
+ Square under tension (quarter model) with mixed Boundary Conditions


### Crack problems in bounded domain
+ Penny-shaped crack under internal pressure in a finite cube 
+ Penny-shaped crack parallel to a free-surface 
+ Bowl-shaped crack // a free-surface  
+ Saw-cut cylindrical sample 


