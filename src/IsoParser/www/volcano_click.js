let interval = 300;
let doubleClickTime = 0;

function singleClickHandler(data) {
  let clickTime = Date.now();

  if ((clickTime - doubleClickTime) > interval) {
    setTimeout(function() {
      if ((clickTime - doubleClickTime) > interval) {
        Shiny.setInputValue("volcano_click", data.points[0].customdata);
      }
    }, interval);
  } else {
    // data is not passed by doubleclick callback
    Shiny.setInputValue("volcano_doubleclick", data.points[0].customdata);
  }
}

function doubleClickHandler() {
  doubleClickTime = Date.now();
}