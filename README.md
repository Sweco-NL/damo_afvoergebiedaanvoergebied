# Generator afwateringseenheden

## Algemeen
Deze python-toolbox is opgezet door Sweco Nederland in opdracht van waterschap Vallei & Veluwe met als doel om uit basisdata afwateringseenheden af te leiden.
Op basis van een raster met een hoogtemodel (maaiveld of grondwaterstand) en een waternetwerk worden stroomgebieden afgeleid, ook wel afwateringseenheden of hydrologische eenheden.
Middels een systematische codering van deze stroomgebieden kan men de stroomgebieden voor verschillende niveaus afleiden.

Er is voor gekozen om open-source te werken. Er wordt geprobeerd om gestructureerd deze python-toolbox op te zetten.

De toolbox bestaat uit een drietal onderdelen die samenwerken, maar ook onafhankelijk zijn te gebruiken:
1. Generator locaties duikers: de overige watergangen (C-watergangen) moeten worden meegenomen, want die zijn (deels) bepalend in de afstroomrichting. Een sluitend netwerk is van belang en zeker op dat niveau zijn de duikers niet bekend (of niet compleet). Middels diverse algoritmes en regels worden de waterdelen en bestaande waternetwerk verbonden. Hierbij wordt de locatie van duikers ingeschat.
2. Generator orde-codering van het netwerk en daarmee voor de afwaterende eenheden (conform [Leidraad Harmoniseren Afvoergebieden](https://kennis.hunzeenaas.nl/file_auth.php/hunzeenaas/a/aa/Leidraden_Harmoniseren_Afvoergebieden_v1.1.pdf)).
3. Generator afwaterende eenheden: op basis van een GHG raster 25x25m de afvoerrichting bepalen en daarmee de afwaterende eenheden. Dit met behulp van o.a. [PyFlwDir van Deltares](https://github.com/Deltares/pyflwdir)

## Installatie environment
We gebruiken pixi om de environment op orde te houden. Installatie van pixi (prefix.dev) kan via de Windows Powershell.
```
iwr -useb https://pixi.sh/install.ps1 | iex
```
Bouw de environment op conform pyproject.toml door vanuit de projectfolder in de command prompt/terminal het volgende te draaien.
```
pixi install
```
We gebruiken ruff voor de code-formatting. Installatie ruff via:
```
pixi global install ruff
```
