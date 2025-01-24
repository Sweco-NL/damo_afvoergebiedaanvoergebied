Voorbeelden
=====================

GeneratorCulvertLocations (workflow Duiker Locaties)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Doel van deze workflow is het verbinden van de C-watergangen met de A/B-watergangen: dit om te bepalen welke gebieden waar aansluiten op de hoofdwatergangen van het waterschap. 
De C-watergangen zijn in beheer van gemeentes of perceeleigenaren en hiervan is vaak geen data beschikbaar. 
Ook vormen de greppels en sloten vaak een verweven netwerk en wil men graag weten welke delen nu waar naartoe afstromen. 

De workflow bestaat uit de volgende stappen:

* Zoekt alle mogelijke verbindingen (duikers) tussen losliggende C-waterdelen en de hoofdwatergangen (A/B);
* Bepaald voor elke mogelijke verbinding/duiker of er een snelweg/spoor/weg of peilgebiedsgrens wordt gekruist;
* Voorziet elke mogelijke verbinding/duiker van een score op basis van criteria (kruisingen, lengte duiker, richting duiker tov watergang);
* Zoekt de duikers met de hoogste scores tussen losliggende C-waterdelen en daarmee ook de verbindingen met de hoofdwatergangen (A/B);
* De richting van de C-watergangen worden gezocht door te bepalen wat de kortste afstand is naar de hoofdwatergangen. Idealiter zou je hier bekijken naar het verhang richting een uitstroompunt naar de hoofdwatergangen.

zie ook `Issue #12 <https://github.com/Sweco-NL/generator_drainage_units/issues/12#issuecomment-2446702722>`_: Selectie beste duiker 

.. image:: _static/generator_culvert_locations_1.jpg
    :alt: Generator Order Levels (workflow duiker-generator)
    :width: 800px
    :align: center

Figuur: duikergenerator - het vinden van de verbindingen van de C-watergangen

.. image:: _static/generator_culvert_locations_2.jpg
    :alt: Generator Order Levels (workflow duiker-generator)
    :width: 800px
    :align: center

Figuur: afleiden stroomrichting C-watergangen naar A/B-watergangen (nu kortste route)


GeneratorOrderLevels (workflow Orde-codering)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Doel van deze workflow is het bepalen van orde nummers en de orde-codering voor iedere watergang en daarmee voor de afwateringseenheden/afvoergebieden die daarmee corresponderen. 
Hiervoor wordt voor de A/B-watergangen uitgegaan van de methode beschreven in de `Leidraad Harmoniseren Afvoergebieden <https://kennis.hunzeenaas.nl/file_auth.php/hunzeenaas/a/aa/Leidraden_Harmoniseren_Afvoergebieden_v1.1.pdf>`_. 

De workflow bestaat (op dit moment) uit de volgende stappen:

* De RWS-wateren waar de watergangen in uitstromen vormen de basis voor de codering (bijv. Veluwemeer: VE, IJssel: IJ);
* De A/B-watergangen die uitstromen in de RWS-wateren zijn van de 2de orde en krijgen een nummer toegewezen die binnen een range ligt die is gespecificeerd voor het waterschap (bijv. Vallei&Veluwe: 712-760). De Leuvenumsebeek krijgt bijvoorbeeld VE.733. Per uitstroompunt zou die vastgelegd moeten worden;
* Ieder individueel watergangsdeel krijgt een opvolgend nummer (bijv. VE.733.001, VE.733.002) of er kan voor gekozen worden dit pas te doen bij splitsingen van de A/B-watergangen;
* Een instromende A/B-watergang wordt een orde hoger (3, 4, 5, etc.) en wordt als gehele zijtak ook meegenomen in de nummering;
* Bij splitsingen of confluences wordt voor het orde nummer (en daarmee ook codering) ervanuit gegaan dat als twee watergangsdelen in het verlengde van elkaar liggen, dat deze van dezelfde orde zijn.
* De C-watergangen die afstromen naar een A/B-watergang worden van een orde hoger dan de watergang waar ze instromen en krijgen dezelfde codering mee (met aanvulling C0001, C0002, ...). Hieruit kan afgeleid worden welke C-watergangen met bijbehorende afvoergebieden bij een watergang horen.

Zie ook: 

* `Issue #16 <https://github.com/Sweco-NL/generator_drainage_units/issues/16#issuecomment-2558479293>`_: Codering RWS wateren en uitstroompunten
* `Issue #17 <https://github.com/Sweco-NL/generator_drainage_units/issues/17#issuecomment-2516835304>`_: Definitie orde A/B watergangen
* `Issue #18 <https://github.com/Sweco-NL/generator_drainage_units/issues/18#issue-2629773652>`_: Definitie orde C watergangen
* `Issue #19 <https://github.com/Sweco-NL/generator_drainage_units/issues/20#issuecomment-2558543651>`_: Definitie orde-codering

.. image:: _static/generator_order_levels_1.jpg
    :alt: Generator Order Levels (workflow orde-codering)
    :width: 800px
    :align: center

Figuur: Afleiden orde nummer van de A/B watergangen

.. image:: _static/generator_order_levels_2.jpg
    :alt: Generator Order Levels (workflow orde-codering)
    :width: 800px
    :align: center

Figuur: Afleiden orde codering van de A/B watergangen

.. image:: _static/order_levels_west_oost.jpg
    :alt: Generator Order Levels (oost)
    :width: 800px
    :align: center

Figuur: Orde nummer van de A/B watergangen voor gehele waterschap Vallei & Veluwe


GeneratorDrainageUnits (workflow Afwateringseenheden)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Workflow voor het genereren van afwateringseenheden: op basis van een GHG raster 25x25m wordt de afvoerrichting bepaald en vervolgens per waterdeel welk gebied erop afstroomt (welke cellen liggen bovenstrooms). 
Er is voor gekozen om te werken met een berekende GHG raster (GHG: gemiddelde hoogste wintergrondwaterstand) omdat dit representatiever is dan het gebruik van het maaiveld, omdat er vooral sprake is van infiltratie en niet oppervlakkige afstroming. 
De analyse gebeurd door middel van een andere open source package `PyFlwDir van Deltares <https://github.com/Deltares/pyflwdir>`_ (Deltares). 
De workflow bestaat (op dit moment) uit de volgende stappen:

* Grof GHG raster wordt gedownscaled naar een opgegeven resolutie.
* Watergangen (lijnen) worden verrasterd. Om te zorgen dat afvoer realistisch richting de watergangen afstroomt wordt het fijne GHG-raster ter hoogte van de watergangen verdiept met 0.20 meter. Deze verlaging wordt minder hoe verder van de watergang.
* Voor het resulterende fijne GHG-raster wordt per cel bepaald welke stroomrichting het heeft (local drainage direction);
* Per watergangsdeel wordt berekend welke cellen er bovenstrooms van liggen. Op de Veluwe kunnen cellen op wel 10-20km afstand liggen die er alsnog naartoe draineren.

.. image:: _static/generator_drainage_units_1.jpg
    :alt: Generator Drainage Units (workflow afwateringseenheden)
    :width: 800px
    :align: center

Figuur: afleiden afwateringseenheden - laaggelegen/polder

.. image:: _static/generator_drainage_units_2.jpg
    :alt: Generator Drainage Units (workflow afwateringseenheden)
    :width: 800px
    :align: center

Figuur: afleiden afwateringseenheden - hogergelegen gebied / vrij-afwaterend

In principe werkt de methode om per watergang het afwaterende gebied te bepalen, alleen de methode kan nog verbeterd worden.
De gebruikte python-package PyFlwDir (net als PCRASTER en vergelijkbare methodes) maakt gebruik van de D8-methode om per cel de afstroomrichting te bepalen aan de hand van de laagste naastliggende cel.

.. image:: _static/ldd_d8.png
    :alt: Generator Drainage Units (ldd d8)
    :width: 300px
    :align: center

Omdat er maar 8 richtingen zijn zie je dit terug in de afwaterende eenheden op aflopende gebieden/hellingen: Bij de Leuvenumsebeek loopt het voornamelijk in de NW richting.
Hier wordt momenteel nog naar gekeken.

Zie ook `Issue #50 <https://github.com/Sweco-NL/generator_drainage_units/issues/50>`_: Omzetten D8-methode naar D-INF-methode.

.. image:: _static/ghg_drainage_units_leuvenumsebeek.jpg
    :alt: Generator Drainage Units (ghg_leuvenumsebeek)
    :width: 800px
    :align: center

Figuur: Leuvenumsebeek, GHG en afwateringseenheden per watergangsdeel


GeneratorNetworkLumping (workflow aggregeren (deel)stroomgebieden)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Workflow om voor opgegeven uitstroompunten het bovenstroomse watersysteem te lumpen en afvoergebieden of (deel)stroomgebieden te genereren. Hierbij wordt overlap gedetecteerd tussen deelstroomgebieden en kan men aangeven hoe de deelgebieden verdeeld worden.
De workflow bestaat (op dit moment) uit de volgende stappen:

* Inladen netwerk van het watersysteem en de bijbehorende afwateringseenheden;
* Definiëren uitstroomlocaties en harde knips in het netwerk;
* Per uitstroompunt zoeken naar gebied bovenstrooms op basis van het netwerk en de richting van de watergangen (deelstroomgebieden);
* Detecteren van overlap tussen deelstroomgebieden en bij welke splitsingen deze gebieden samen komen;
* Voor deze splitsingen bepalen welke richting prioriteit heeft;
* Deelstroomgebieden afronden door afwateringseenheden eraan te koppelen.


