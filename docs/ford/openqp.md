---
project: OPENQP
version: {!VERSION!}
project_github: https://github.com/Open-Quantum-Platform/openqp
project_download: https://github.com/Open-Quantum-Platform/openqp/releases/tag/release
summary: ![OPENQP](./media/logo/openqp-logo.png)
         {: style="text-align: center" }
author: OpenQP Team
author_description:
github: https://github.com/Open-Quantum-Platform/
email: cheolho.choi@gmail.com
docmark: +
extra_filetypes: pro ;
                 py  #
preprocess: false
media_dir: @CMAKE_SOURCE_DIR@/docs/ford/media
md_base_dir: @CMAKE_SOURCE_DIR@
page_dir: @CMAKE_SOURCE_DIR@/docs/ford/guide
src_dir: @CMAKE_SOURCE_DIR@/source
output_dir: @CMAKE_SOURCE_DIR@/doc/ford
md_extensions: markdown.extensions.toc
source: true
graph: true
search: false
favicon: @CMAKE_SOURCE_DIR@/docs/ford/media/logo/logo.png

---

Open Quantum Platform (OpenQP) is a quantum chemical platform featuring cutting-edge capabilities
     like Mixed-Reference Spin-Flip (MRSF)-TDDFT, with an emphasis on an open-source ecosystem.

For more information, check out the [user guide](https://github.com/Open-Quantum-Platform/openqp/wiki).

This is the documentation for the main branch of OpenQP.

## Citing OpenQP
If you use OpenQP in your research, please cite the following papers:

```
@article{mironov2024openqp,
  title={OpenQP: A Quantum Chemical Platform Featuring MRSF-TDDFT with an Emphasis on Open-Source Ecosystem},
  author={Mironov, Vladimir and Komarov, Konstantin and Li, Jingbai and Gerasimov, Igor and Nakata, Hiroya and Mazaherifar, Mohsen and Ishimura, Kazuya and Park, Woojin and Lashkaripour, Alireza and Oh, Minseok and others},    journal={Journal of Chemical Theory and Computation},
  volume={20},
  number={21},
  pages={9464--9477},
  year={2024},
  publisher={ACS Publications}
}
```

```
@article{park2023mixed,
  title={Mixed-reference spin-flip time-dependent density functional theory: Multireference advantages with the practicality of linear response theory},
  author={Park, Woojin and Komarov, Konstantin and Lee, Seunghoon and Choi, Cheol Ho},
  journal={The Journal of Physical Chemistry Letters},
  volume={14},
  number={39},
  pages={8896--8908},
  year={2023},
  publisher={ACS Publications}
}
```

```
@article{lee2018eliminating,
  title={Eliminating spin-contamination of spin-flip time dependent density functional theory within linear response formalism by the use of zeroth-order mixed-reference (MR) reduced density matrix},
  author={Lee, Seunghoon and Filatov, Michael and Lee, Sangyoub and Choi, Cheol Ho},
  journal={The Journal of chemical physics},
  volume={149},
  number={10},
  year={2018},
  publisher={AIP Publishing}
}
```

```
@article{lee2019efficient,
  title={Efficient implementations of analytic energy gradient for mixed-reference spin-flip time-dependent density functional theory (MRSF-TDDFT)},
  author={Lee, Seunghoon and Kim, Emma Eunji and Nakata, Hiroya and Lee, Sangyoub and Choi, Cheol Ho},
  journal={The Journal of Chemical Physics},
  volume={150},
  number={18},
  year={2019},
  publisher={AIP Publishing}
}
```
