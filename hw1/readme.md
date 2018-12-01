Protein ligand binding affinity is one of the many methods available for probing a protein&#39;s binding affinity. By mutating a protein and measuring its binding affinity to its native ligand, a scientist can deduce how different alterations to a protein's sequence can alter protein structure, and deploy that information to engineer proteins with higher/lower binding affinities to a target molecule.

In a simple protein ligand system, using mass action kinetics we recover the following system:

[E]+[S]k<sub>1</sub>⇌k<sub>-1</sub>[ES]

where E represents the enzyme, S the substrate, and ES the enzyme substrate complex. The formation constant is represented by k<sub>1</sub>, and disassociation k<sub>-1</sub>.

From this single equation, we can derive the following ODE using mass action kinetics:
d[ES]dt=k<sub>1</sub>[E][S]−k<sub>-1</sub>[ES]

Assuming steady state:

0=k<sub>1</sub>[E][S]−k<sub>-1</sub>[ES]

[ES][E][S]=k<sub>1</sub>k<sub>-1</sub>=KD−1

Where K<sub>D</sub> is a general disassociation constant. In the whole system, since there are [E]<sub>t</sub> total enzymes in the system, and enzymes and substrates bind in a one to one ratio, we can write the following:

[E]<sub>t</sub>=[E]+[ES] 

where [E] represents the concentration of unassociated enzyme in the system.

The fraction of bound enzyme to total enzyme in the system is

[ES][E]<sub>t</sub>=[E][S]K<sub>D</sub>[E]+[E][S]K<sub>D</sub>=[S]K<sub>D</sub>+[S]

Which results in the equation we'll use to model our system.


