"use strict";

/*----------------------------------------------------------------------------*/

function initialize() {
    if (window.plugin) {
        plugin.destroy();
    }
    LiteMolWrapper();
    window.plugin = LiteMol.Plugin.create({target: "#litemol", viewportBackground: "#FFFFFF", layoutState: {hideControls: true}});
}

initialize();

/*----------------------------------------------------------------------------*/

document.querySelector("#loadDomainButton").onclick = function() {
    var inputDomainFile = document.querySelector("#domainInput").files[0];
    var reader = new FileReader();
    reader.onload = drawDomain.bind(this, reader);
    reader.readAsText(inputDomainFile);
}

function drawDomain(reader) {
    var domain = JSON.parse(reader.result);
    var pdbId = domain.pdbId;
    var colorMap = new Map([
        ["skipped", "skipped"],
        ["d", {"r": 0.839, "g": 0.910, "b": 0.976}],
        ["m", {"r": 0.749, "g": 0.937, "b": 0.561}],
        ["i", {"r": 0.749, "g": 0.937, "b": 0.561}],
        ["0", {"r": 0.946392, "g": 0.930761, "b": 0.442367}],
        ["1", {"r": 0.983196, "g": 0.743758, "b": 0.138453}],
        ["2", {"r": 0.977092, "g": 0.55085, "b": 0.03905}],
        ["3", {"r": 0.916462, "g": 0.387481, "b": 0.164924}],
        ["4", {"r": 0.801871, "g": 0.258674, "b": 0.283099}],
        ["5", {"r": 0.658463, "g": 0.178962, "b": 0.372748}],
        ["6", {"r": 0.497257, "g": 0.119379, "b": 0.424488}],
        ["7", {"r": 0.3415, "g": 0.062325, "b": 0.429425}],
        ["8", {"r": 0.169575, "g": 0.042489, "b": 0.340874}],
        ["9", {"r": 0.029432, "g": 0.021503, "b": 0.114621}]]);
    var colorMask = domain.colorMask.map(function (i) {return colorMap.get(i);});
    plugin.destroy();
    LiteMolWrapper(colorMask);
    window.plugin = LiteMol.Plugin.create({target: "#litemol", viewportBackground: "#FFFFFF", layoutState: {hideControls: true}});
    plugin.loadMolecule({id: pdbId, url: "https://www.ebi.ac.uk/pdbe/static/entry/" + pdbId + ".cif", format: "cif"});
}
