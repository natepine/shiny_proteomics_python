let correlationInterval = 300;
let correlationDoubleClickTime = 0;

function correlationSingleClickHandler(data) {
  let clickTime = Date.now();

  if ((clickTime - correlationDoubleClickTime) > correlationInterval) {
    setTimeout(function() {
      if ((clickTime - correlationDoubleClickTime) > correlationInterval) {
        Shiny.setInputValue("correlation_click", data.points[0].customdata);
      }
    }, correlationInterval);
  } else {
    // data is not passed by doubleclick callback
    Shiny.setInputValue("correlation_doubleclick", data.points[0].customdata);
  }
}

function correlationDoubleClickHandler() {
  correlationDoubleClickTime = Date.now();
}