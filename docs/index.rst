.. DrainageUnits documentation master file, created by
   sphinx-quickstart on Wed Dec 18 09:19:01 2024.


Generator Drainage Units: Workflows voor hydrologische basisdata waterschappen
============================================================================================================================
o.a. workflows voor voorbewerking waternetwerk, verbinden c-watergangen middels duikers, afleiden van afwateringseenheden en stroomgebieden, automatische toekennen van orde-codering. 


Algemeen
----------------------------
Deze python-toolbox is opgezet door Sweco Nederland binnen twee losse opdrachten voor waterschap Vallei & Veluwe en Aa & Maas met als doel om op basis van hydrologische basisdata van het watersysteem analyses uit te voeren. Hierbij wordt open-source gewerkt aan deze public repository waarbij we ons richten op overzichtelijkheid en voldoende documentatie en voorbeeld-scripts inclusief testdata.


Waterschap Vallei & Veluwe
----------------------------
De vraag om op basis van een raster met een hoogtemodel (maaiveld of in dit geval grondwaterstand GHG) en gegeven waternetwerk stroomgebiedjes af te leiden, ook wel afwateringseenheden of hydrologische eenheden genoemd. Middels netwerk-analyse en codering is het mogelijk om deze te aggregeren tot elk gewenst niveau:

- **generator_culvert_locations**: workflow die voortbouwt op een al bestaande 'duikergenerator' van het waterschap waarin de locaties van duikers voor de C-watergangen worden bepaald. Dit gebeurd op basis van (configureerbare) regels, die rekening houden met kruizingen van (spoor)wegen en peilgebiedsgrenzen, de lengte van de duiker (hoe korter, hoe beter) en de richting van de duiker ten opzichte van de watergang (zelfde hoek heeft voorkeur). 

- **generator_order_levels**: workflow voor het bepalen van orde nummers en de orde-codering van iedere watergang en daarmee voor de afwaterende eenheden (conform `Leidraad Harmoniseren Afvoergebieden <https://kennis.hunzeenaas.nl/file_auth.php/hunzeenaas/a/aa/Leidraden_Harmoniseren_Afvoergebieden_v1.1.pdf>`_), op basis van de codering kan eenvoudig geaggregeerd worden. Ook worden hier de C-watergangen (niet hoofdwatergangen) meegenomen in de analyse.

- **generator_drainage_units**: workflow voor het genereren van afwateringseenheden: op basis van een GHG raster 25x25m de afvoerrichting bepalen en daarmee de afwaterende eenheden. Dit met behulp van andere open source packages zoals `PyFlwDir van Deltares <https://github.com/Deltares/pyflwdir>`_


Waterschap Aa & Maas
----------------------------
De vraag om op basis van benedenstroomse uitstroompunten (deel)stroomgebieden te genereren.

- **generator_network_lumping**: workflow om voor opgegeven uitstroompunten het bovenstroomse watersysteem te lumpen en afvoergebieden of (deel)stroomgebieden te genereren. Hierbij wordt overlap gedetecteerd en kan men aangeven hoe de deelgebieden verdeeld worden.


.. toctree::
   :maxdepth: 2
   :caption: Inhoud:
   
   Installatie <getting_started>
   Voorbeelden <some_examples>
   API docs <api_docs>



