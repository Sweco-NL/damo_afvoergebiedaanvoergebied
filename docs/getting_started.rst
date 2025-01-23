Getting Started
=====================

Voor gebruik download de public repository. Daarin zit zowel de package generator_drainage_units als enkele voorbeeldscripts inclusief testdata.


Installatie environment
----------------------------
We gebruiken `pixi <https://pixi.sh/>`_ (prefix.dev) om de environment op orde te houden. Installatie van pixi kan via de Windows Powershell::

    iwr -useb https://pixi.sh/install.ps1 | iex

Bouw de environment op op basis van het bestand pyproject.toml door in de projectfolder via de terminal te draaien::

    pixi install

Vervolgens wordt een pixi environment opgebouwd in de folder .pixi/envs (waarschijnlijk default-environment).


Klaarzetten basisdata
----------------------------
Er moet een locatie van een map worden opgegeven waar gewerkt wordt. Hierbinnen komt een map met basisdata en een map met resultaten. 

De basisdata bestaat (voor nu) uit geopackages met steeds een layer/laag met dezelfde naam als het geopackage. Het raster bestand voor de GHG is van het format TIF/ASC/NetCDF:
Bij het aanroepen van de analyse wordt de benodigde data uit de map met basisdata ingelezen en als dit vereist is ook uit de resultaten-folder (resultaten van duikergenerator wordt gebruikt voor orde-codering en vervolgens voor het genereren van de afwaterende eenheden)
De basis-data bestaat uit:

Datasets voor het watersysteem:

- **rws_wateren.gpkg**: Polygonen die de RWS-wateren voorstellen inclusief code
- **hydroobjecten.gpkg**: A/B-watergangen (hoofdwatergangen waterschappen)
- **hydroobjecten_extra.gpkg**: Extra hoofdwatergangen optioneel toe te voegen
- **overige_watergangen.gpkg**: C-watergangen (lijn-elementen afgeleid van waterdelen)

Datasets die gebruikt worden in de analyse van de 'Duikergenerator' (om te bepalen of er lijn-elementen gekruist worden):

- **keringen.gpkg**: waterkeringen
- **nwb.gpkg**: centerlines van alle wegen (nationaal wegen bestand)
- **snelwegen.gpkg**: centerlines van de snelwegen
- **spoorwegen.gpkg**: centerlines van de spoorwegen
- **peilgebieden.gpkg**: peilgebieden (polygonen) om de peilgebeidsgrenzen mee te nemen.

Topografische hoogtemodel (raster data) voor afleiden afwateringseenheden. In dit geval GHG-raster, maar kan ook maaiveld (DTM):

- **GHG_2000-2010_L1.NC** Het rasterbestand kan apart worden opgegeven en ingeladen.

Afleiden deelstroomgebieden middels uitstroompunten:

- **inflow_outflow_points.gpkg**: punten die dienen als uitstroompunt (eventueel ook als inlaatpunten te gebruiken).


Draaien van de voorbeeld notebooks
----------------------------
Voor het draaien van de notebooks, gebruiken wij standaard VS Code. Als kernel kan de python-executable in de map .pixi/envs/default/python.exe worden geselecteerd. Middels Ctrl+Shift+P en "Python: Select Interpreter" kan de kernel geselecteerd worden.

De analyse inclusief brondata en de resultaten komt in een map terecht. In de notebooks wordt hierbij verwezen naar de .env-file. Dit is gedaan omdat een ieder de data op een ander path kan hebben staan.

