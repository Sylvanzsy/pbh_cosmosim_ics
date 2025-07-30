# PHANTOM: Primordial black Holes And Nonlinear perTurbations fOr siMulations

This document describes the Fortran code `pbh_ic`, designed to insert Primordial Black Holes (PBHs) into cosmological initial conditions (ICs) formatted for Gadget simulations. The code supports various customizable options to modify particle distributions, perturbations, and streaming velocities relevant for cosmological simulations.

## Code Structure and Functionality

The main goal of the program is to read a pre-existing Gadget-format IC file, insert PBHs into the particle data, apply relevant perturbations, and output the modified IC file.

### Configuration Flags

The following preprocessor flags control code functionality:

* `DEBUG`: Enables verbose debugging output for development purposes.
* `OUTPBH`: Writes an additional output file containing initial positions and velocities of PBHs.
* `PERTURB`: Enables perturbations of particles by PBHs.
* `STREAMING`: Adds relative streaming velocities between baryons and dark matter (DM).
* `FLIPID`: Swaps particle IDs between gas and dark matter particles (useful for specific codes like Arepo).
* `PBH`: Activates insertion of PBHs into the simulation.
* `DAMPING`: Introduces damping to perturbations.
* `VELOCITY`: Includes velocity perturbations.
* `TRUNCATE`: Applies truncation to particle displacements.
* `CENTER_PBH`: Places a single PBH at the simulation box center.
* `ALLPBH`: Perturbs particles using all PBHs rather than nearest ones only.
* `PERIOD`: Uses periodic boundary conditions for calculating distances.
* `NONLINEAR`: Introduces initial nonlinear structures around PBHs.
* `PERTURB_PBH`: Enables self-perturbations among PBHs.
* `NO_GAS_PERTURB`: Disables perturbations to gas particles.
* `noPBH`: Includes only perturbations induced by PBHs without explicitly inserting PBHs.

### Important Parameters

The code defines several physical and numerical parameters crucial for PBH insertion:

* `mpbh`: PBH mass (in solar masses, \$M\_\odot\$).
* `fpbh`: Mass fraction of DM in PBHs.
* `sigma_rms`: Root mean square streaming velocity at recombination (default 0.8).
* `z_rec`: Redshift of recombination (default 1100).
* `a`: Scale factor corresponding to the IC redshift.

## Physics and Cosmology

### Growth of Perturbations

The displacement (`\psi`) and velocity perturbations induced by PBHs are computed based on cosmological linear and non-linear perturbation theory. Here, we adapt the analytical form from ([Jiao et al. (2023)](https://ui.adsabs.harvard.edu/abs/2023PhRvD.108d3510J)) for perturbations induced by point objects. Key equations include:

#### Linear displacement

$$
\psi_{l}(r, a) = \frac{3}{2}\frac{GM_{\mathrm{PBH}}t_0^2}{r^2}\left[0.4\left(\frac{a_{\mathrm{eq}}}{a}\right)^{3/2} + 0.6\frac{a}{a_{\mathrm{eq}}} - 1\right] \frac{f_{\mathrm{scale}}}{\mathrm{kpc/h}}
$$

#### Linear velocity perturbation

$$
\dot{\psi_{l}}(r, a) = \frac{0.6 GM_{\mathrm{PBH}}t_0}{a_{\mathrm{eq}}^{3/2}r^2}\left[\left(\frac{a_{\mathrm{eq}}}{a}\right)^{1/2} - \left(\frac{a_{\mathrm{eq}}}{a}\right)^3\right]\frac{a f_{\mathrm{scale}}}{10^5 \mathrm{km/s}}
$$

### Non-linear Turnaround Radius

The non-linear scale (`q_nl`) for turnaround radius around a PBH is given by:

$$
q_{nl}(M_{\mathrm{PBH}}, a) = \left(\frac{1.8 GM_{\mathrm{PBH}} t_0^2 a}{a_{\mathrm{eq}}}\right)^{1/3}
$$

### Non-linear displacement and velocity

For particles within the turnaround radius:

$$
\psi_{nl}(r,a) = r\left[1-\frac{a_{ta}(r)}{4a}\right], \quad \dot{\psi}_{nl}(r,a)=\frac{2}{3}\frac{\psi(r)}{t},
$$

where \$a\_{ta}(r)\$ is the turnaround scale factor defined by:

$$
a_{ta}(r)=a_{eq}\frac{5r^3}{9GM_{\mathrm{PBH}}t_0^2}
$$

For detailed theoretical background, see relevant cosmological structure formation literature such as ([Press & Schechter (1974)](https://ui.adsabs.harvard.edu/abs/1974ApJ...187..425P)) and subsequent PBH formation papers like e.g. ([Inman & Ali-Haïmoud (2019)](https://ui.adsabs.harvard.edu/abs/2019PhRvD.100h3528I)).

## Algorithmic Workflow

1. **Read IC File:** Reads existing particle positions, velocities, IDs, and masses.
2. **Streaming Velocity:** Adds relative velocity between gas and DM particles if enabled.
3. **PBH Placement:** Inserts PBHs either randomly or at the center of the box.
4. **Perturbation Calculation:** Computes gravitational perturbations due to PBHs.
5. **Particle Displacement:** Modifies positions and velocities based on linear/non-linear perturbation theory.
6. **Output:** Writes new IC file including PBHs.

## Output File Naming Convention

The output IC file name indicates selected options:

* `_pbh`: PBHs included.
* `_pbhcen`: PBH placed at the center.
* `_streaming`: DM-baryon streaming velocity included.
* `_ngasptb`: No perturbation to gas.

## Example Usage

To compile and run the code using the provided Makefile:

```bash
make clean
make
./pbh_ic.exe
```

Ensure the correct Gadget IC file is referenced in the `path` and `icbase` parameters.

## Dependencies

* `random_object`: Random number generation module.
* `kdtree`: Efficient nearest neighbor search.

## Citation

When using this code in scientific work, please cite the following work:

```
@ARTICLE{2022MNRAS.514.2376L,
       author = {{Liu}, Boyuan and {Zhang}, Saiyang and {Bromm}, Volker},
        title = "{Effects of stellar-mass primordial black holes on first star formation}",
      journal = {\mnras},
     keywords = {black hole physics, dark ages, reionization, first stars, dark matter, early Universe, Astrophysics - Astrophysics of Galaxies, Astrophysics - Cosmology and Nongalactic Astrophysics},
         year = 2022,
        month = aug,
       volume = {514},
       number = {2},
        pages = {2376-2396},
          doi = {10.1093/mnras/stac1472},
archivePrefix = {arXiv},
       eprint = {2204.06330},
 primaryClass = {astro-ph.GA},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2022MNRAS.514.2376L},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}



@ARTICLE{2023arXiv231204085L,
       author = {{Liu}, Boyuan and {Bromm}, Volker},
        title = "{Impact of primordial black holes on the formation of the first stars and galaxies}",
      journal = {arXiv e-prints},
     keywords = {Astrophysics - Astrophysics of Galaxies, Astrophysics - Cosmology and Nongalactic Astrophysics},
         year = 2023,
        month = dec,
          eid = {arXiv:2312.04085},
        pages = {arXiv:2312.04085},
          doi = {10.48550/arXiv.2312.04085},
archivePrefix = {arXiv},
       eprint = {2312.04085},
 primaryClass = {astro-ph.GA},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2023arXiv231204085L},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}


@ARTICLE{2024ApJ...975..139Z,
       author = {{Zhang}, Saiyang and {Bromm}, Volker and {Liu}, Boyuan},
        title = "{How Do Primordial Black Holes Change the Halo Mass Function and Structure?}",
      journal = {\apj},
     keywords = {Cosmological perturbation theory, Early universe, Dark matter, Large-scale structure of the universe, 341, 435, 353, 902, Astrophysics - Cosmology and Nongalactic Astrophysics, Astrophysics - Astrophysics of Galaxies},
         year = 2024,
        month = nov,
       volume = {975},
       number = {1},
          eid = {139},
        pages = {139},
          doi = {10.3847/1538-4357/ad7b0d},
archivePrefix = {arXiv},
       eprint = {2405.11381},
 primaryClass = {astro-ph.CO},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2024ApJ...975..139Z},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```



For theoretical background, also cite relevant cosmological structure formation literature, e.g., ([Inman & Ali-Haïmoud (2019)](https://ui.adsabs.harvard.edu/abs/2019PhRvD.100h3528I)), and Gadget code for reading the initial conditions: ([Hopkins & Raives (2015)](https://ui.adsabs.harvard.edu/abs/2016MNRAS.455...51H)).

This tool provides a robust method to study cosmological impacts of PBHs through numerical simulations, aiding research into early structure formation and cosmological constraints on dark matter models.
