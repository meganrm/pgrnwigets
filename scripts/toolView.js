
var $grid = $('.grid').isotope({
  itemSelector: '.element-item',

});


// filter items on button click
$('.filter-button-group').on('click', '.btn-filter', function() {
  $('.about').empty();

  var filterValue = $(this).attr('data-filter');
  $grid.isotope({ filter: filterValue });
  var cat = $(this).attr('data-filter').replace(/\./g, '')
  var catInfo = Tool.allCategories.filter(function(ele){
    return ele.categoryID === cat;
  })[0]
  console.log(cat, catInfo);
  var source = $('#category-template').html();
  var renderTemplate = Handlebars.compile(source);
  $('.about').append(renderTemplate(catInfo));
});

// quick search regex
var qsRegex;

// use value of search field to filter
var $quicksearch = $('.quicksearch').keyup(debounce(function() {
  qsRegex = new RegExp($quicksearch.val(), 'gi' );
  $grid.isotope({  filter: function() {
      return qsRegex ? $(this).attr("class").match( qsRegex ) : true;
    }});
}, 200));

// debounce so filtering doesn't happen every millisecond
function debounce( fn, threshold ) {
  var timeout;
  return function debounced() {
    if ( timeout ) {
      clearTimeout( timeout );
    }
    function delayed() {
      fn();
      timeout = null;
    }
    timeout = setTimeout( delayed, threshold || 100 );
  }
}

$('[data-filter=".topChoice"]').click();
