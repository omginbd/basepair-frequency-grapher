/*******************
 *
 * file: stringified version of uploaded file
 *
 * returns json representation of reads:
 * {
 *  "onemers": [
 *    {
 *      "A":
 *        {
 *          "count": 12345,
 *          "proportion": 0.75
 *        },
 *      ...
 *    },
 *    ...
 *  ],
 *  "twomers": [...],
 *  ...
 * }
 *
 *******************/
const parseFile = file => {
  file = file.replace(/(Position: [0-9]*)|(GC content: [0-9]*.[0-9]*)/g, "");
  var lines = file.split("\n");
  lines = lines.map((line) => line.trim());
  lines = lines.map((line) => line.split(';'));

  var parsedObj = new Object();
  parsedObj.onemers   = new Array();
  parsedObj.twomers   = new Array();
  parsedObj.threemers = new Array();
  parsedObj.fourmers  = new Array();

  var unparsedOnemers   = new Array();
  var unparsedTwomers   = new Array();
  var unparsedThreemers = new Array();
  var unparsedFourmers  = new Array();

  lines.forEach(function(line) {
    var currLine = line[0];
    if (currLine.startsWith("A=")) {
      unparsedOnemers.push(line);
    } else if (currLine.startsWith("AA=")) {
      unparsedTwomers.push(line);
    } else if (currLine.startsWith("AAA=")) {
      unparsedThreemers.push(line);
    } else if (currLine.startsWith("AAAA=")) {
      unparsedFourmers.push(line);
    }
  });

  // Parse onemers
  unparsedOnemers.forEach(function(onemer) {
    // make a onemer object
    var newOnemer = new Object();
    // for every grouping in onemer
    onemer.forEach(function(grouping, index, arr) {
      var parsedGrouping = grouping.replace(/(=\[)|(:)/g, ";").replace(/(\])/g, "").split(";");
      newOnemer[parsedGrouping[0]] = new Object();
      newOnemer[parsedGrouping[0]] = {
        count: Number(parsedGrouping[1]),
        proportion: Number(parsedGrouping[2])
      };
    });
    // put the new onemer in the parsed object
    parsedObj.onemers.push(newOnemer);
  });

  // Parse twomers
  unparsedTwomers.forEach(function(twomer) {
    // make a twomer object
    var newTwomer = new Object();
    // for every grouping in twomer
    twomer.forEach(function(grouping, index, arr) {
      var parsedGrouping = grouping.replace(/(=\[)|(:)/g, ";").replace(/(\])/g, "").split(";");
      newTwomer[parsedGrouping[0]] = new Object();
      newTwomer[parsedGrouping[0]] = {
        count: Number(parsedGrouping[1]),
        proportion: Number(parsedGrouping[2])
      };
    });
    // put the new twomer in the parsed object
    parsedObj.twomers.push(newTwomer);
  });

  // Parse threemers
  unparsedThreemers.forEach(function(threemer) {
    // make a threemer object
    var newThreemer = new Object();
    // for every grouping in threemer
    threemer.forEach(function(grouping, index, arr) {
      var parsedGrouping = grouping.replace(/(=\[)|(:)/g, ";").replace(/(\])/g, "").split(";");
      newThreemer[parsedGrouping[0]] = new Object();
      newThreemer[parsedGrouping[0]] = {
        count: Number(parsedGrouping[1]),
        proportion: Number(parsedGrouping[2])
      };
    });
    // put the new threemer in the parsed object
    parsedObj.threemers.push(newThreemer);
  });

  // Parse fourmers
  unparsedFourmers.forEach(function(fourmer) {
    // make a fourmer object
    var newFourmer = new Object();
    // for every grouping in fourmer
    fourmer.forEach(function(grouping, index, arr) {
      var parsedGrouping = grouping.replace(/(=\[)|(:)/g, ";").replace(/(\])/g, "").split(";");
      newFourmer[parsedGrouping[0]] = new Object();
      newFourmer[parsedGrouping[0]] = {
        count: Number(parsedGrouping[1]),
        proportion: Number(parsedGrouping[2])
      };
    });
    // put the new fourmer in the parsed object
    parsedObj.fourmers.push(newFourmer);
  });

  return parsedObj;
}

/*******************
 *
 * read1: json representation of file1
 * read2: json representation of file2
 * (see above for data format)
 *
 * returns multidimensional array of p-values:
 * [
 *  [.1, .2, .3, ...],
 *  [...],
 *  ...
 * ]
 *
 *******************/
const compareReads = (read1, read2) => {}

const beginProcessing = () => {
  const file1Container = document.getElementById("file1")
  const file2Container = document.getElementById("file2")
  const file1Reader = new FileReader()
  const file2Reader = new FileReader()
  file1Reader.readAsText(file1Container.files[0])
  file2Reader.readAsText(file2Container.files[0])
  awaitProcessing(file1Reader, file2Reader)
}

const awaitProcessing = (f1, f2) => {
  if (f1.readyState === 2 && f2.readyState === 2) {
    updateTables(compareReads(parseFile(f1.result), parseFile(f2.result)))
  } else {
    setTimeout(awaitProcessing, 100, f1, f2)
  }
}

const updateTables = pValueTable => {
  for (let key in pValueTable) {
    let tableString = "<table><tbody>"
    for (let j = 0; j < pValueTable[key].length; j++) {
      tableString += "<tr>"
      for (let k = 0; k < pValueTable[key][j]; k++) {
        tableString += `<td>${pValueTable[key][j][k]}</td>`
      }
      tableString += "</tr>"
    }
    tableString += "</tbody></table>"
    const renderNode = document.getElementById(`table${key}`)
    renderNode.innerHTML = tableString
  }
}

// Init function
;(function() {
  const goButton = document.getElementById("go")
  goButton.addEventListener("click", beginProcessing)
})()
