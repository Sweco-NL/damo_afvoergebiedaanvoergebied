# Generator afwateringseenheden

Deze python-toolbox is opgezet door Sweco Nederland met als doel om uit basisdata afwateringseenheden af te leiden.
Op basis van een raster met een hoogtemodel (maaiveld of grondwaterstand) en een waternetwerk worden stroomgebieden afgeleid, ook wel afwateringseenheden of hydrologische eenheden.
Middels een systematische codering van deze stroomgebieden kan men de stroomgebieden voor verschillende niveaus afleiden.

Dit wordt uitgevoerd in opdracht van waterschap Vallei & Veluwe. Er is voor gekozen om open-source te werken. 
Er wordt daarom ook geprobeerd om gestructureerd deze python toolbox op te zetten.

De toolbox bestaat uit een drietal onderdelen die samenwerken, maar ook onafhankelijk zijn te gebruiken:
1. Generator locaties duikers: we willen de overige watergangen (C-watergangen) meenemen, deze zijn deels bepalend in de afstroomrichting. Op dat niveau zijn echter de duikers niet bekend, alleen de waterdelen. Middels diverse algoritmes en regels wordt een verbinding gezocht middels de meest logische duikers.
2. Generator afwaterende eenheden: op basis van GHG raster 25x25m willen we op basis van afvoerrichting de afwaterende eenheden afleiden. Dit middels de python-moduel [PyFlwDir van Deltares](https://github.com/Deltares/pyflwdir)
3. Generator orde-codering van het netwerk of afwaterende eenheden (conform [Leidraad Harmoniseren Afvoergebieden](https://kennis.hunzeenaas.nl/file_auth.php/hunzeenaas/a/aa/Leidraden_Harmoniseren_Afvoergebieden_v1.1.pdf)).
