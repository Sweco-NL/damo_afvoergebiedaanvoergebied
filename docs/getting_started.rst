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

Draaien van de voorbeeld notebooks
----------------------------
Voor het draaien van de notebooks, gebruiken wij standaard VS Code. Als kernel kan de python-executable in de map .pixi/envs/default/python.exe worden geselecteerd. Middels Ctrl+Shift+P en "Python: Select Interpreter" kan de kernel geselecteerd worden.

