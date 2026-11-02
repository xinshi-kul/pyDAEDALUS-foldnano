# pyDAEDALUS for Foldnano

A [Foldnano](https://github.com/XinShiLab/foldnano) plugin that routes a wireframe mesh as scaffolded DNA origami with two-helix (DX) edges. It is adapted from [pyDAEDALUS](https://github.com/lcbb/pyDAEDALUS) by the Bathe BioNanoLab, MIT, with the authors' permission, under the original GNU General Public License, version 2 or any later version.

The plugin is not listed in Foldnano's plugin Hub yet. Wireframe regions are a developer preview in Foldnano, available while plugin dev mode is on.

## What it does

Run it on an empty wireframe region with a reference mesh. It routes the mesh with DAEDALUS (B-form DNA, M13 scaffold, the minimum edge length you choose) and adds the helices, scaffold and staples to the region. DAEDALUS sets the size of the structure from the edge length, so the plugin rescales the region's mesh to match.

The mesh must be closed: every edge shared by exactly two faces. `PLY_Files/` holds the upstream example meshes, which all meet that rule; open one in Foldnano to try the plugin.

## Branches

- `master` tracks the upstream project unchanged.
- `foldnano_plugin` holds the plugin. The routing code in `pyDAEDALUS/Automated_Design/` is upstream's, with one NumPy 2 fix. The server, web interface, PDB output and example scripts are not part of this branch. `pyDAEDALUS/NOTICE` lists every change.

## Layout

`pyDAEDALUS/` holds exactly the files of the plugin zip: `fplugin.json` (the plugin manifest), `plugin.py` (the Foldnano entry point), `convert_cndo.py` (turns the routing into Foldnano helices and strands), `Automated_Design/` (the routing code), `LICENSE.txt` and `NOTICE`.

## Releasing a version

1. Change the plugin files and bump `version` in `pyDAEDALUS/fplugin.json`.
2. Pack the zip with the [foldnano-plugins](https://github.com/XinShiLab/foldnano-plugins) repo, checked out next to this one:

   ```sh
   cd ../foldnano-plugins
   npm run pack -- ../pyDAEDALUS-foldnano/pyDAEDALUS
   ```

   Packing is reproducible: the same files always give the same SHA-256, so anyone can check a listed zip against this branch.

3. Add the printed entry to the listing there, following its `CONTRIBUTING.md`.

## How to cite

If you use this plugin in published work, please cite:

```
Veneziano, R., Ratanalert, S., Zhang, K., Zhang, F., Yan, H., Chiu, W. and Bathe, M.
Designer nanoscale DNA assemblies programmed from the top down.
Science 352(6293), 1534 (2016).
https://doi.org/10.1126/science.aaf4388
```

## License

GNU General Public License, version 2 or any later version, as the original project. See `LICENSE.txt`.
