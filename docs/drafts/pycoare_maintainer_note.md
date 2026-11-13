# Draft note to the pycoare maintainer

*A draft for Chris to edit and send (as a GitHub issue on pyCOARE/coare, or by email).
Not sent.*

**Subject: A sea-ice extension built on pycoare — courtesy note and two questions**

Hello,

I'm a researcher at the University of East Anglia working on air-sea-ice fluxes in the
Greenland Sea. With Ian Renfrew I'm building a small open-source package,
`coare-seaice`, that extends COARE to sea ice and the marginal ice zone: form drag after
Lüpkes et al. (2012) and Elvidge et al. (2016), the Andreas (1987) scalar roughness over
ice, and an ice surface temperature from a surface energy balance. It uses pycoare for the
open-water fluxes and depends only on its public interface.

Repository: https://github.com/cbarrell16/coare-seaice (MIT, early development).

I wanted to let you know before releasing anything, and to ask two things:

1. **Name.** Are you comfortable with `coare-seaice`? I chose it to avoid looking like an
   official part of pycoare, and I'm happy to change it.
2. **Stable interface.** We read the final roughness length, the 10 m neutral wind, the
   gust factor and the friction velocity from the `coare_35` / `coare_36` objects. Are
   those attributes, and the functions in `pycoare.util`, something you consider stable
   across 0.x releases?

The modules are written so they could be contributed to pycoare later if that were of
interest to you. Thank you for maintaining pycoare.

Chris Barrell
