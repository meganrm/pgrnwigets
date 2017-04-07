var ref = new Firebase("https://pgrnqualtricsdate.firebaseio.com");
var userresponse = ref.child('responses');

var fromfirebase=[]

userresponse.once("value", function(snapshot) {
  snapshot.forEach(function(childSnapshot) {
    var key = childSnapshot.key();
    var childData = childSnapshot.val();
    if (childSnapshot.child('status').val()=='approved'){
      fromfirebase.push(childData)
    }
  });
  renderMetaData('#resources', createMetaData('resources'))
  renderMetaData('#expertise', createMetaData('expertise'));
  renderMetaData('#dieases', createMetaData('dieases'));
  renderMetaData('#research-state', researchState);


});


function newTagjQuery(tagType, parent, text, newclass){
  var elementName=$('<'+tagType+'>');
  elementName.appendTo(parent).text(text).addClass(newclass);
  return elementName;
}

function newlistitem(listid, text, id) {
  var li= newTagjQuery('li', listid, text);
  li.attr(id);
  return li;
}


userresponse.on('child_added', function (snapshot) {
    var results = snapshot.val();
    if (results['status']=='approved'){
		var firstname = results['firstname'];
		var lastname = results['lastname'];
	    var name= firstname +' '+lastname;
	    fullname=results['fullname']
		// var li= newlistitem(membersList, ' ',   fullname);
		// newTagjQuery('span', li,  fullname, 'name')
	  //   newTagjQuery('span', li,  results['display_intitution'], 'affiliation')
	  }

})

function renderMetaData($listId, array){
  var source = $('#hist-template').html();
  var renderTemplate = Handlebars.compile(source);
  for (var i = 0; i < array.length; i++) {
    if (array[i].count > 1 ) {
      $(renderTemplate(array[i])).appendTo($($listId))
    }
  }
}

function createMetaData(cat) {
  allResources = fromfirebase.map(function(ele){return ele[cat]}).filter(function(ele){return ele.length>0})
  hist = allResources.reduce(function(acc, cur){cur.split(',').forEach(function(li){
  acc.push(li.trim())}); return acc},[])
  .reduce(function(acc2, cur2) {
      var repeat = false
      for (var i = 0; i < acc2.length; i++) {
        if (acc2[i].name === cur2) {
          repeat = i
          break
        }
      }
    if (parseInt(repeat)) {
      acc2[repeat].count ++;
    } else {
      obj = {};
      obj.name = cur2;
      obj.count = 1;
      acc2.push(obj);
      }
      return acc2
    },[])

  return hist.sort(function(a, b){return b.count - a.count})
}

//TODO: print txt
