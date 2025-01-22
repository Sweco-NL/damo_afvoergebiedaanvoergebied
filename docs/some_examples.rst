Voorbeelden
=====================

GeneratorCulvertLocations (workflow Duiker Locaties)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Doel van deze workflow is het verbinden van de C-watergangen met de A/B-watergangen: om te bepalen welke gebieden stromen in de hoofdwatergangen van het waterschap, dat is namelijk vaak niet duidelijk of het is verweven.
De workflow bestaat uit de volgende onderdelen:

* zoekt alle mogelijke verbindingen (duikers) tussen losliggende C-waterdelen en de hoofdwatergangen (A/B);
* bepaald voor elke mogelijke verbinding/duiker of er een snelweg/spoor/weg of peilgebiedsgrens wordt gekruisd;
* voorziet elke mogelijke verbinding/duiker van een score op basis van criteria (kruisingen, lengte duiker, richting duiker tov watergang);
* zoekt de duikers met de hoogste scores tussen losliggende C-waterdelen en daarmee ook de verbindingen met de hoofdwatergangen (A/B);
* de richting van de C-watergangen worden gezocht door te bepalen wat de kortste afstand is naar de hoofdwatergangen.

.. image:: assets/generator_culvert_locations.png
    :alt: Generator Culvert Locations (generator duiker locaties)
    :width: 600px
    :align: center


GeneratorOrderLevels (workflow Orde-codering)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Generates order levels and order codes for all hydroobjects including lower level hydroobjects

* punt 1
* punt 2
* punt 3

.. image:: assets/generator_order_levels.png
    :alt: Generator Order Levels (generator orde-codering)
    :width: 600px
    :align: center


GeneratorDrainageUnits (workflow Orde-codering)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Generates drainage units for each hydroobject based on a terrain model

* punt 1
* punt 2
* punt 3

.. image:: assets/generator_drainage_units.png
    :alt: Generator Drainage Units (generator afwateringseenheden)
    :width: 600px
    :align: center


GeneratorNetworkLumping (workflow genereren (deel)stroomgebieden)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Generates upstream (sub)basins for predefined outflow points

* punt 1
* punt 2
* punt 3

.. image:: assets/generator_network_lumping.png
    :alt: Generator Network Lumping (generator stroomgebieden)
    :width: 600px
    :align: center

