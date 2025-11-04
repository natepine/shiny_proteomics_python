Shiny.addCustomMessageHandler('clearTooltips', function(empty) {
   $('.tooltip').remove();
});

Shiny.addCustomMessageHandler('resetValue', function(variableName) {
   Shiny.onInputChange(variableName, null);
});

Shiny.addCustomMessageHandler('setLabelText', function(msg) {
   const [id, lbl] = msg;
   $("label[for='" + id + "']").text(lbl);
});

