function Tool (opts) {
  for (var key in opts) {
    if (opts[key] !== 'NA') {
      this[key] = opts[key];
    }
  }
  if (opts) {
    this.keywordsList = opts.Keywords.split(', ');
    this.categoryID = opts.Category.replace(/\s/g, '');
  }
}

Tool.allTools = []

Tool.allCategories = tools
  .map(function(ele){
    return {
      categoryID : ele.Category.replace(/\s/g, ''),
      Category : ele.Category,
      subtitle : ele.Subtitle
    }
  }).filter(function(ele, index, array){
    return array.map(function(mapItem){ return mapItem['categoryID']; }).indexOf(ele['categoryID']) === index;
})

Tool.renderAll = function(templateid, parent, array) {
  var source = $(templateid).html();
  var renderTemplate = Handlebars.compile(source);
  for (var i = 0; i < array.length; i++) {
    $(parent).append(renderTemplate(array[i]))
  }
}

Tool.loadAll = function loadAll(array) {
  var googlekeys = ['name', 'url', 'Keywords', 'Category', 'Subtitle', 'stats', 'Descriptions']
  var encodedArray = []
  for (var j = 1; j < array.length; j++) {
    var row = array[j]
    var rowObj = {}
    for (var k = 0; k < row.length; k++) {
      rowObj[googlekeys[k]] = row[k]
    }
    tool = new Tool(rowObj)
    encodedArray.push(tool)
  }
  return encodedArray
}

Tool.fetchAllGoogle = function () {
  url = 'https://sheets.googleapis.com/v4/spreadsheets/1VXFAgVBJ3NU0TdWS0GFWh1VYbpcuQEvFgaBEh9nqFeA/values/Sheet1!A:G?key=AIzaSyBw6HZ7Y4J1dATyC4-_mKmt3u0hLRRqthQ';
  return new Promise(function (resolve, reject) {
    $.get(url, function (response) {
      var rows = response.values
      if (rows.length === 0) {
        console.error('No data found')
      } else {
        var encodedArray = Tool.loadAll(rows)
        console.log(encodedArray);
        resolve(encodedArray)
      }
    })
  })
}


$(window).resize(function(){
  var footerHeight = $('.footer').height();
  var headerHeight = $('.banner-wrap').height();
  var navHeight = $('.birdseye-header').height();
  var docHeight = window.innerHeight;
  $('iframe').height(docHeight - navHeight - headerHeight - footerHeight)
})

function loadTools() {
  Tool.fetchAllGoogle().then(function(tools){
    Tool.renderAll('#toolGrid-template', '.grid', tools)
    Tool.renderAll('#buttons-template', '.dropdown-menu', Tool.allCategories)
    $grid = $('.grid').isotope({
      itemSelector: '.element-item',
    });
  })
}

loadTools()
