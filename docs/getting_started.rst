Installatie en gebruik van de workflows
=====================

Download voor gebruik de public repository. Daarin zit zowel de package damo_afvoergebiedaanvoergebied als enkele voorbeeldscripts inclusief testdata.


Installatie environment
----------------------------
We gebruiken `pixi <https://pixi.sh/>`_ (prefix.dev) om de environment op orde te houden. Installatie van pixi kan via de Windows Powershell::

    iwr -useb https://pixi.sh/install.ps1 | iex

Bouw de environment op op basis van het bestand pyproject.toml door in de projectfolder via de terminal te draaien::

    pixi install

Vervolgens wordt een pixi environment opgebouwd in de folder .pixi/envs (waarschijnlijk default-environment).


Klaarzetten basisdata
----------------------------
Er wordt gewerkt vanuit een vaste map die opgegeven dient te worden bij iedere analyse. Hierbinnen komt een map met basisdata en een map met resultaten. De naam van deze sub-folders kan opgegeven worden. 

De basisdata hoort te bestaan uit geopackages met steeds een layer/laag met dezelfde naam als het geopackage. Het raster bestand voor de GHG is van het format TIF/ASC/NetCDF.
Bij het aanroepen van de analyse wordt de benodigde data uit de map met basisdata ingelezen en als dit vereist is ook uit de resultaten-folder (resultaten van duikergenerator wordt gebruikt voor orde-codering en vervolgens voor het genereren van de afwaterende eenheden)
De basis-data bestaat uit:

Datasets voor het watersysteem:

- **rws_water.gpkg**: Polygonen van RWS-wateren inclusief officiële code
- **hydroobject.gpkg**: A-/B-watergangen (hoofdwatergangen waterschappen)
- **hydroobject_extra.gpkg**: Extra hoofdwatergangen optioneel toe te voegen
- **overige_watergang.gpkg**: C-watergangen (lijn-elementen afgeleid van waterdelen)

Datasets die gebruikt worden in de analyse van de 'Duikergenerator' (om te bepalen of er lijn-elementen gekruist worden):

- **kering.gpkg**: waterkering
- **weg.gpkg**: centerlines van alle wegen (nationaal wegen bestand)
- **snelweg.gpkg**: centerlines van de snelweg
- **spoorweg.gpkg**: centerlines van de spoorweg
- **peilgebied.gpkg**: peilgebied (polygonen) om de peilgebeidsgrenzen mee te nemen.

Topografische hoogtemodel (raster data) voor afleiden afwateringseenheden. In dit geval GHG-raster, maar kan ook maaiveld (DTM):

- **GHG.NC** Het rasterbestand kan apart worden opgegeven en ingeladen.

Afleiden deelstroomgebieden middels uitstroompunten:

- **outflow_nodes_hydro.gpkg**: punten die worden opgegeven als uitstroompunt waarvoor dan stroomgebieden worden afgeleid.
- **splitsing.gpkg**: Voor splitsingen in het netwerk kan worden opgegeven welke richting voorrang heeft, zodat overlap voorkomen wordt. Deze locaties worden bij elke analyse ook gedetecteerd en weggeschreven.


Gebruik van de (voorbeeld-)notebooks
----------------------------
Voor het draaien van de notebooks, gebruiken wij standaard VS Code. Als kernel kan de python-executable in de map .pixi/envs/default/python.exe worden geselecteerd (zie boven). Middels Ctrl+Shift+P en "Python: Select Interpreter" kan de kernel geselecteerd worden.

De analyse inclusief brondata en de resultaten komt in een map terecht. In de notebooks wordt hierbij verwezen naar de .env-file. Dit is gedaan omdat een ieder de data op een ander path kan hebben staan.

De verschillende workflows kunnen los gedraaid worden, maar er zit ook een volgorde in:

- Workflow Duikers: bereid het watersysteem helemaal voor, voert voorbewerking uit op de A-/B-watergangen, verbindt de C-watergangen middels duikers, etc; 
- Workflow GebiedsOrde: Genereert order-nummers en orde-codering voor alle (A-/B-/C-)watergangen; 
- Workflow Afvoergebieden: Genereert de afwateringseenheden (polygonen); 
- Workflow NetworkLumping: Genereert de deelstroomgebieden voor uitstroompunten, (sluitend) waternetwerk en afwateringseenheden.

Om de analyses uit te voeren zijn voor iedere analyse losse functies geschreven:

- `run_generator_duikers() <https://generator-drainage-units.readthedocs.io/en/latest/api_docs.html#damo_afvoergebiedaanvoergebied.generator_duikers.run_generator_duikers>`_
- `run_generator_gebiedsorde() <https://generator-drainage-units.readthedocs.io/en/latest/api_docs.html#damo_afvoergebiedaanvoergebied.generator_gebiedsorde.run_generator_gebiedsorde>`_
- `run_damo_afvoergebiedaanvoergebied() <https://generator-drainage-units.readthedocs.io/en/latest/api_docs.html#damo_afvoergebiedaanvoergebied.generator_network_lumping.run_generator_network_lumping>`_
- `run_generator_network_lumping() <https://generator-drainage-units.readthedocs.io/en/latest/api_docs.html#damo_afvoergebiedaanvoergebied.generator_network_lumping.run_generator_network_lumping>`_

Alle functies hebben eigenlijk alleen het path nodig waar de map met basisdata staat. De resultaatmap wordt automatisch aangemaakt. Indien gewenst kunnen extra parameters worden meegegeven voor het draaien van de functies.
