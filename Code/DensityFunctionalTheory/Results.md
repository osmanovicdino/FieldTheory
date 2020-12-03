We wish to demonstrate that the MIPs system we have defined leads to Phase separation. 

For a system with the pair potential:

$$ H(r)=\exp(-r^{10}/(2d)) $$

Where r is the distance between two particles. These particles have a non-conservative force in $(x,y)$ of

$$F_{nc} =(\cos(z),\sin(z))$$

For any pair of particles, this gives us an effective potential for the particle (a smeared out version of this)

We wish to test whether a system with this free energy displays phase separation. To that level we construct an effective free energy for this system:

$$ F[\rho] = \int \mathrm{d}r f_{hs}(\rho) $$

Where 

$$f_{hs}(\rho)=\rho\big(\ln(\rho \lambda^3)-1+F^{HS}[\rho_{2D}]\big)$$

where the hard disk term depends only on the 2D density

In scaled particle theory, the compressibility of disks is given by:

$$\frac{1}{(1-y)^2}$$

Using the condition at the bottom, the excess free energy goes as:

$$ F^{HS}= \frac{y}{1-y}-\log(1-y) $$

where y is the two dimensional density: $y=\int\mathrm{d}z \rho (x,y,z)\pi d^2/4$

where we take the functional derivative with respect to $\rho$

The functional derivative of the ideal term is given by:
$$ \log(\rho)$$

the effective free energy due to the interactions we have introduced are given by:

$$F_{int}= \iint \mathrm{d}r\mathrm{d}r'\phi(r,r')\rho(r)\rho(r')$$

So the derivative of this term is given by:

$$ \int \phi(r,r')\rho(r')\mathrm{d}r'$$

This leads to the minimization condition:

$$\log(\rho)+\frac{\log (1-y(\rho ))-y(\rho )}{(1-y(\rho ))^2}+\int \phi(r,r')\rho(r')\mathrm{d}r'=0$$
## Links

### Hard Disk Theory:

http://www.sklogwiki.org/SklogWiki/index.php/Scaled-particle_theory#Equation_of_state_of_hard_disks

### Relationship between free energy and compressibility:

from https://sci-hub.se/https://aip.scitation.org/doi/10.1063/1.4997256


$$ \frac{\partial}{\partial \rho}\bigg(\frac{\beta F}{N}\bigg)=\frac{Z}{\rho}$$


