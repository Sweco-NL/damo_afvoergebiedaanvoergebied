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
