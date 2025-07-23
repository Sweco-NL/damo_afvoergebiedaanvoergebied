## Workflows network-analyses op hydrologische basisdata waterschappen

o.a. voorbewerking waternetwerk, verbinden c-watergangen middels duikers, afwateringseenheden, stroomgebieden afleiden, automatische orde-codering. Bekijk de uitgebreide documentatie via de [read-the-docs](https://generator-drainage-units.readthedocs.io/en/latest/)


### Installatie environment
We gebruiken pixi om de environment op orde te houden. Installatie van pixi (prefix.dev) kan via de Windows Powershell:
```
iwr -useb https://pixi.sh/install.ps1 | iex
```
Bouw de environment op op basis van het bestand pyproject.toml door in de projectfolder via de terminal te draaien:
```
pixi install
```

### Voorbereiding project
Om de verschillende generators te kunnen draaien, moet er in de lokale clone van de repository een bestand `.env` worden toegevoegd:

```
📁 generator_drainage_units
├── ..
└── .env
```

Deze bevat een `BASE_DIR`, de locatie waar per `case_name` data `basisdata` wordt opgeslagen en `resultaat` wordt gegenereerd. De inhoud van `.env` wordt daarmee bijvoorbeeld:
```
BASE_DIR="d:\pad\naar\mijn\analyses"
```

De binnen de folder `BASE_DIR` zit een sub-folder `case_name`, bijvoorbeeld `Leuvenumse_beek`, waaronder in een map `0_basisdata` alle gegevens zijn opgeslagen die nodig zijn voor het uitvoeren van de generatoren. In `1_resultaat` komen alle resultaten.

```
📁 BASE_DIR
└── 📁 case_name
    ├── 📁 0_basisdata
    │   ├── GHG_2000-2010_L1.NC
    │   ├── hydroobjecten.gpkg
    │   ├── inflow_outflow_points.gpkg
    │   ├── inflow_outflow_splits.gpkg
    │   ├── keringen.gpkg
    │   ├── nwb.gpkg
    │   ├── overige_watergangen.gpkg
    │   ├── peilgebieden.gpkg
    │   ├── rws_wateren.gpkg
    │   ├── snelwegen.gpkg
    │   ├── spoorwegen.gpkg
    └── 📁 1_resultaat
        ├── ...
```

### Genereren uitvoer

In het mapje `notebooks` vindt je drie notebooks waarmee je de generatoren kunt draaien. Voor een beschrijving van de generatoren zie de [documentatie](https://generator-drainage-units.readthedocs.io/en/latest/description_workflows.html#generatorculvertlocations-workflow-duiker-locaties).

```
📁 generator_drainage_units
├── ...
├── 📁 notebooks
│   ├── run_generator_culvert_locations.ipynb
│   ├── run_generator_order_code.ipynb
│   └── run_generator_drainage_units_d16.ipynb
├── ...
```
