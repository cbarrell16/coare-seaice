# Draft note to the pycoare maintainer

*A draft for Chris to edit and send. Not sent.*

*Maintainer (from pycoare's `pyproject.toml`): Andrew Scherer, Oregon State University.
Send either as a GitHub issue on pyCOARE/coare (public) or by email to the address listed
there. Post the bug report (`pycoare_issue_rh_mutation.md`) first, so this note can link
to it. Send this before the name is registered on PyPI.*

**Subject: A sea-ice extension built on pycoare — courtesy note and two questions**

Hello Andrew,

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

While building on pycoare we found that `qair` changes the caller's `rh` array in
place; I've opened a separate issue with a reproduction (link).

The modules are written so they could be contributed to pycoare later if that were of
interest to you. Thank you for maintaining pycoare.

Chris Barrell
