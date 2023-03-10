# FlutterBridge_FD
The coupled flutter velocity of a single-span suspension bridge is computed in the frequency domain



## Summary
The critical flutter velocity Vcr of a suspension bridge is estimated using a simple computational model accounting for the lateral, vertical and torsional motion of the bridge deck and using a multimodal approach. The computation is conducted in the frequency domain using the method proposed in [1]. The computed value of Vcr is compared to the famous analytical expressions from Selberg [2] and Rocard [3].


## Content

The present submission contains:

- The function flutterFd, which computes the critical flutter velocity following [1]
- The function VcrFlutter, which computes the critical flutter velocity following [2,3]
- An example file Example.m
- Two .mat files modalParameters_case1.mat and modalParameters_case2.mat that are used to load the eigenfrequencies and mode shapes of the two bridge models investigated

This is the first version of the submission. Some bugs may still exist. Any question, comments or suggestion is warmly welcomed.

## References:

[1] Jain, A., Jones, N. P., & Scanlan, R. H. (1996). Coupled aeroelastic and aerodynamic response analysis of long-span bridges. Journal of Wind Engineering and Industrial Aerodynamics, 60, 69-80.

[2] Selberg, A., & Hansen, E. H. (1966). Aerodynamic stability and related aspects of suspension bridges.

[3] Rocard, Y. (1963). Instabilite des ponts suspendus dans le vent-experiences sur modele reduit. Nat. Phys. Lab. Paper, 10.
