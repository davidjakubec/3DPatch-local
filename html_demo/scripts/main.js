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
    var colorMask = domain.colorMask;
    plugin.destroy();
    LiteMolWrapper(colorMask);
    window.plugin = LiteMol.Plugin.create({target: "#litemol", viewportBackground: "#FFFFFF", layoutState: {hideControls: true}});
    plugin.loadMolecule({id: pdbId, url: "https://www.ebi.ac.uk/pdbe/static/entry/" + pdbId + "_updated.cif", format: "cif"});
}
