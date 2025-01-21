Voorbeelden
=====================

GeneratorCulvertLocations (workflow Duiker Locaties)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Doel van deze module is het verbinden van de C-watergangen met de A/B-watergangen: om te bepalen welke gebieden stromen in de hoofdwatergangen van het waterschap, dat is namelijk vaak niet duidelijk of het is verweven.
De module bestaat uit de volgende onderdelen:
- zoekt alle mogelijke verbindingen (duikers) tussen losliggende C-waterdelen en de hoofdwatergangen (A/B);
- bepaald voor elke mogelijke verbinding/duiker of er een snelweg/spoor/weg of peilgebiedsgrens wordt gekruisd;
- voorziet elke mogelijke verbinding/duiker van een score op basis van criteria;
- zoekt de meest logische verbindingen (met de hoogste score) tussen losliggende C-waterdelen en daarmee ook de verbindingen met de hoofdwatergangen (A/B);
- de richting van de C-watergangen worden gezocht door te bepalen wat de kortste afstand is naar de hoofdwatergangen.

GeneratorOrderLevels (workflow Orde-codering)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Generates order levels and order codes for all hydroobjects including lower level hydroobjects

GeneratorDrainageUnits (workflow Orde-codering)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Generates drainage units for each hydroobject based on a terrain model

GeneratorNetworkLumping (workflow genereren (deel)stroomgebieden)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Generates upstream (sub)basins for predefined outflow points

Utils (help-functies)
^^^^
Several help-functions