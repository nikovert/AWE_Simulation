# AWE_Simulation

This Repo provides both the code used for controller synthesis as well as the simulation framework used for testing. The accompanying paper is available as a preprint  [here](https://www.vertovec.info/publication/safety-aware-hybrid-control-of-airborne-wind-energy-systems/)

## Abstract

In order to move Airborne Wind Energy (AWE) into the future, a fundamental concern is being able to guarantee that safety requirements placed on the systems are met.  Due to the high dimensional complexity of AWE systems, however, strict mathematical robustness guarantees become difficult to compute. 

We draw on research from high dimensional Hamilton-Jacobi (HJ) reachability analysis to compute the optimal trajectory for tracking a figure-eight flight path during the pumping cycle of AWE systems while enforcing safety constraints on the system, such as those placed on the aircraft speeds and the tether force. In addition to providing the optimal control policy, the subzero level-set of the computed value function inherent in HJ reachability analysis indicates the backward reachable set (BRS), the set of states from which it is possible to safely drive the system into a target set within a given time without entering undesirable states, defined by an avoid set.

The BRS can then be used as a maneuverability envelope to derive a switching law, such that the safety controller can be used in conjunction with arbitrary least restrictive controllers to provide a safe hybrid control law. In such a setup, the safety controller is only needed when the system approaches the boundary of the maneuverability envelope. Such a hybrid control law is a notable improvement over existing robust control approaches that assume the worst-case environmental and system behavior at all times, leading to potentially sub-optimal control laws.

## Getting Started

### Dependencies
The Simulink model requires an installation of Simulink as well as a couple of Matlab toolboxes, including the Aerospace Toolbox. For controller synthesis, we use some tools from the helperOC toolbox which is an extension of the level-set methods toolbox by Ian Mitchell. Only the level-set methods toolbox needs to be installed, which can be found here: ('https://www.cs.ubc.ca/~mitchell/ToolboxLS/'). 

### Executing program

To generate the control lookup table used by Simulink, you need to run the following script
```
hj_reachability/awe/main.m
```
Depending on your computer specs, you made need to adjust the grid size to account for memory constraints. Note that a too small grid, will result in an inaccurate BRS calculation, which might make the controller underperform. During execution, a file called tables.mat will be created, which will be used by the Simulink Simulation.

Having generated the control lookup table, the Simulation can be run by executing 
```
simulink/run_sim.m
```

## Authors

The coding for this piece of research was done by Nikolaus Vertovec, and builds on top of the work of Berkley Hybrid Systems Lab, which developed the level-set-methods toolbox as well as the doctoral work of Sebastian Rapp (http://awesco.eu/project/esr02/).

## License

This project is licensed under the GNU License - see the LICENSE.md file for details

## Acknowledgments

Inspiration, code snippets, etc.
* [helperOC](https://github.com/HJReachability/helperOC)
* [TU Delft Repo](https://data.4tu.nl/articles/software/Scripts_for_AWE_Control_Design_and_Simulation/13172666)
