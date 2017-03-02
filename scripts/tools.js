function Tool (opts) {
  for (var key in opts) {
    this[key] = opts[key];
    this.keywordsList = opts.Keywords.split(', ');
    this.categoryID = opts.Category.replace(/\s/g, '');
  }
}

Tool.allTools = []

Tool.allCategories = tools
  .map(function(ele){return ele.Category.replace(/\s/g, '')}).filter(function(ele, index, array){
  return array.indexOf(ele) === index;
})

// init Isotope


Tool.renderAll = function(templateid, parent, array) {
  var source = $(templateid).html();
  var renderTemplate = Handlebars.compile(source);
  for (var i = 0; i < array.length; i++) {
    if (typeof(array[i]) !== 'string') {
      item = new Tool(array[i])
      Tool.allTools.push(item)
    } else {
      item = array[i]
    }
    $(parent).append(renderTemplate(item))
  }
}

Tool.renderAll('#toolGrid-template', '.grid', tools)
Tool.renderAll('#buttons-template', '.dropdown-menu', Tool.allCategories)
