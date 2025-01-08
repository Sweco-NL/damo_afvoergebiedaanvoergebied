.. DrainageUnits documentation master file, created by
   sphinx-quickstart on Wed Dec 18 09:19:01 2024.


Generator Drainage Units: Workflows network-analyses op hydrologische basisdata waterschappen
============================================================================================================================
o.a. voorbewerking waternetwerk, verbinden c-watergangen middels duikers, afwateringseenheden, stroomgebieden afleiden, automatische orde-codering. Bekijk de uitgebreide documentatie via de `read-the-docs <https://generator-drainage-units.readthedocs.io/>`


Algemeen
----------------------------
Deze python-toolbox is opgezet door Sweco Nederland binnen twee losse opdrachten voor waterschap Aa en Maas en Vallei & Veluwe met als doel om uit hydrologische basisdata van de waterschappen netwerk-analyses uit te voeren. Er is voor gekozen om open-source te werken. We proberen deze python-toolbox op te zetten inclusief testdata, voorbeeld-scripts en documentatie.


Waterschap Vallei & Veluwe
----------------------------
De vraag om op basis van een raster met een hoogtemodel (maaiveld of grondwaterstand) en een waternetwerk stroomgebiedjes af te leiden, ook wel afwateringseenheden of hydrologische eenheden genoemd. Middels codering zou het mogelijk moeten zijn om deze te aggregeren tot elk gewenst niveau:

- **generator_culvert_locations**: voortbouwend op een al bestaande 'duikergenerator' worden de locaties van duikers bepaald. Dit gebeurd op basis van (configureerbare) regels, die rekening houden met kruizingen van wegen en peilgebiedsgrenzen, de lengte van de duiker (hoe lager, hoe beter) en de richting van de duiker ten opzichte van de watergang (zelfde hoek heeft voorkeur). 

- **generator_order_levels**: bepalen van de orde en de orde-codering van het netwerk en daarmee voor de afwaterende eenheden (conform `Leidraad Harmoniseren Afvoergebieden <https://kennis.hunzeenaas.nl/file_auth.php/hunzeenaas/a/aa/Leidraden_Harmoniseren_Afvoergebieden_v1.1.pdf>`)

- **generator_drainage_units**: workflow voor het genereren van afwateringseenheden: op basis van een GHG raster 25x25m de afvoerrichting bepalen en daarmee de afwaterende eenheden. Dit met behulp van o.a. `RESPHIGI <https://gitlab.com/deltares/imod/respighi>` en `PyFlwDir van Deltares <https://github.com/Deltares/pyflwdir>`


Waterschap Aa en Maas
----------------------------
De vraag om op basis van benedenstroomse uitstroompunten (deel)stroomgebieden te genereren.

- **generator_network_lumping**: Toolbox om voor gegeven uitstroompunten het bovenstroomse netwerk te lumpen en afvoergebieden of (deel)stroomgebieden te genereren.


Installatie environment
----------------------------
We gebruiken pixi om de environment op orde te houden. Installatie van pixi (prefix.dev) kan via de Windows Powershell:

```
iwr -useb https://pixi.sh/install.ps1 | iex
```
Bouw de environment op op basis van het bestand pyproject.toml door in de projectfolder via de terminal te draaien:

```
pixi install
```

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
   API docs <preprocess>



