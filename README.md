# ellipseCollisionDetection
Determines if 2 ellipses collide/overlap using a simple approach.

The final mathematical criteria are used from Etayo et al., 2006;
A new approach to characterizing the relative position of two ellipses depending on one parameter:
http://www.sciencedirect.com/science/article/pii/S0167839606000033?via%3Dihub

Two programs exist in this repository:
> a ROOT CERN C++ file
* root .L overlapEllipses.cxx
* overlapEllipses(angleA, semiMajorA, semiMinorA, hA, kA, angleB, semiMajorB, semiMinorB, hB, kB, 0 or 1)
* A refers to ellipse A, B to ellipse B, the last function parameter will display some information if set to 1
> a Mathematica notebook
* Open the notebook in Wolfram Mathematica
* Change the parameters at the top
* Evaluation (toolbar) -> Evaluate Notebook

Each program determines if 2 ellipses overlap, or not, and then sketches the results.

Example output:
![alt text](https://github.com/LukeBatten/ellipseCollisionDetection/blob/master/img/ellipseOverlap1.png)
![alt text](https://github.com/LukeBatten/ellipseCollisionDetection/blob/master/img/ellipseOverlap2.png)
