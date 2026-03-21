## Workflow voor hydrologische basisdata waterschappen

o.a. voorbewerking waternetwerk, verbinden c-watergangen middels duikers, afwateringseenheden, automatische orde-codering, stroomgebieden afleiden. Bekijk de uitgebreide documentatie via de [read-the-docs](https://damo-afvoergebiedaanvoergebied.readthedocs.io/nl/latest/)


DOCUMENTATIE WORDT IN MAART NOG AANGEPAST!
AFGELOPEN TIJD IS ER VEEL GEWERKT OM ALLERLEI KLEINE BUGS ERUIT TE HALEN.
DAARNAAST HEBBEN WE DE OVERSTAP GEMAAKT NAAR EEN TOML-BESTAND ALS INPUT WAAR ALLE SETTINGS IN STAAN. DAT MAAKT HET VEEL OVERZICHTELIJKER, MAAR HIERDOOR LOOPT DE DOCUMENTATIE ACHTER.


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
Om de verschillende workflows te kunnen draaien, moet er in de lokale clone van de repository een bestand `.env` worden toegevoegd:

```
📁 damo_afvoergebiedaanvoergebied
├── ..
└── .env
```

Deze bevat een `BASE_DIR`, de locatie waar per `case_name` data `basisdata` wordt opgeslagen en `resultaat` wordt gegenereerd. De inhoud van `.env` wordt daarmee bijvoorbeeld:
```
BASE_DIR="d:\pad\naar\mijn\analyses"
```

De binnen de folder `BASE_DIR` zit een sub-folder `case_name`, bijvoorbeeld `Leuvenumse_beek`, waaronder in een map `0_basisdata` alle gegevens zijn opgeslagen die nodig zijn voor het uitvoeren van de workflows. In `1_resultaat` komen alle resultaten.

```
📁 BASE_DIR
└── 📁 case_name
    ├── 📁 0_basisdata
    │   ├── ghg.nc
    │   ├── hydroobject.gpkg
    │   ├── kering.gpkg
    │   ├── weg.gpkg
    │   ├── overige_watergang.gpkg
    │   ├── peilgebied.gpkg
    │   ├── rws_water.gpkg
    │   ├── snelweg.gpkg
    │   ├── spoorweg.gpkg
    └── 📁 1_resultaat
        ├── ...
```

### Genereren uitvoer

In het mapje `notebooks` vindt je drie notebooks waarmee je de workflows kunt draaien. Voor een beschrijving van de workflows zie de [documentatie](https://damo-afvoergebiedaanvoergebied.readthedocs.io/en/latest/description_workflows.html#GeneratorDuikers-workflow-duiker-locaties).

```
📁 damo_afvoergebiedaanvoergebied
├── ...
├── 📁 notebooks
│   ├── run_generator_duikers.ipynb
│   ├── run_generator_order_code.ipynb
│   └── run_generator_afvoergebieden_d16.ipynb
├── ...
```
