Shiny.addCustomMessageHandler('resetValue', function(variableName) {
   Shiny.onInputChange(variableName, null);
});

Shiny.addCustomMessageHandler('setDatasetName', function(name) {
   $("body div.wrapper header.main-header nav.navbar").append("<div id=datasetName class=ellipsis>" + name + "</div>")
});

Shiny.addCustomMessageHandler('setPTM', function(empty){
   $('.skin-blue .main-header .logo, .skin-blue .main-header .navbar').css('background-color', '#58BD6A');
   $('.skin-blue .main-header .logo, .skin-blue .main-header .navbar .sidebar-toggle').hover(
      function(){
         $(this).css('background-color', '#333232');
      },
      function(){
         $(this).css('background-color', '#58BD6A');
      }
   );
   $('.skin-blue .main-header .logo').text('TMT Mosaic PTM');
   $('a[href$=\"#shiny-tab-summary\"] > span').text('Site Summary');
});

