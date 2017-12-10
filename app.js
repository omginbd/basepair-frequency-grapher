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
  file = file.replace(/(Position: [0-9]*)|(GC content: [0-9]*.[0-9]*)/g, "")
  var lines = file.split("\n")
  lines = lines.map(line => line.trim())
  lines = lines.map(line => line.split(";"))

  var parsedObj = new Object()
  parsedObj.onemers = new Array()
  parsedObj.twomers = new Array()
  parsedObj.threemers = new Array()
  parsedObj.fourmers = new Array()

  var unparsedOnemers = new Array()
  var unparsedTwomers = new Array()
  var unparsedThreemers = new Array()
  var unparsedFourmers = new Array()

  lines.forEach(function(line) {
    var currLine = line[0]
    if (currLine.startsWith("A=")) {
      unparsedOnemers.push(line)
    } else if (currLine.startsWith("AA=")) {
      unparsedTwomers.push(line)
    } else if (currLine.startsWith("AAA=")) {
      unparsedThreemers.push(line)
    } else if (currLine.startsWith("AAAA=")) {
      unparsedFourmers.push(line)
    }
  })

  // Parse onemers
  unparsedOnemers.forEach(function(onemer) {
    // make a onemer object
    var newOnemer = new Object()
    // for every grouping in onemer
    onemer.forEach(function(grouping, index, arr) {
      var parsedGrouping = grouping
        .replace(/(=\[)|(:)/g, ";")
        .replace(/(\])/g, "")
        .split(";")
      newOnemer[parsedGrouping[0]] = new Object()
      newOnemer[parsedGrouping[0]] = {
        count: Number(parsedGrouping[1]),
        proportion: Number(parsedGrouping[2])
      }
    })
    // put the new onemer in the parsed object
    parsedObj.onemers.push(newOnemer)
  })

  // Parse twomers
  unparsedTwomers.forEach(function(twomer) {
    // make a twomer object
    var newTwomer = new Object()
    // for every grouping in twomer
    twomer.forEach(function(grouping, index, arr) {
      var parsedGrouping = grouping
        .replace(/(=\[)|(:)/g, ";")
        .replace(/(\])/g, "")
        .split(";")
      newTwomer[parsedGrouping[0]] = new Object()
      newTwomer[parsedGrouping[0]] = {
        count: Number(parsedGrouping[1]),
        proportion: Number(parsedGrouping[2])
      }
    })
    // put the new twomer in the parsed object
    parsedObj.twomers.push(newTwomer)
  })

  // Parse threemers
  unparsedThreemers.forEach(function(threemer) {
    // make a threemer object
    var newThreemer = new Object()
    // for every grouping in threemer
    threemer.forEach(function(grouping, index, arr) {
      var parsedGrouping = grouping
        .replace(/(=\[)|(:)/g, ";")
        .replace(/(\])/g, "")
        .split(";")
      newThreemer[parsedGrouping[0]] = new Object()
      newThreemer[parsedGrouping[0]] = {
        count: Number(parsedGrouping[1]),
        proportion: Number(parsedGrouping[2])
      }
    })
    // put the new threemer in the parsed object
    parsedObj.threemers.push(newThreemer)
  })

  // Parse fourmers
  unparsedFourmers.forEach(function(fourmer) {
    // make a fourmer object
    var newFourmer = new Object()
    // for every grouping in fourmer
    fourmer.forEach(function(grouping, index, arr) {
      var parsedGrouping = grouping
        .replace(/(=\[)|(:)/g, ";")
        .replace(/(\])/g, "")
        .split(";")
      newFourmer[parsedGrouping[0]] = new Object()
      newFourmer[parsedGrouping[0]] = {
        count: Number(parsedGrouping[1]),
        proportion: Number(parsedGrouping[2])
      }
    })
    // put the new fourmer in the parsed object
    parsedObj.fourmers.push(newFourmer)
  })

  return parsedObj
}

/*******************
 *
 * read1: json representation of file1
 * read2: json representation of file2
 * (see above for data format)
 *
 * returns multidimensional array of p-values:
 * {
 *  onemers: [
 *    {
 *      A: pvalue,
 *      T,
 *      ...
 *    },
 *    ...
 *  ],
 *  twomers: [],
 *  ...
 * }
 *
 *******************/
const compareReads = (parsed1, parsed2) => {
  // z = ((p1 - p2) - 0) / (Math.sqrt(po * (1 - po) * (1/n1 + 1/n2)))
  var pValueObjects = new Object()
  pValueObjects.onemers = new Array()
  pValueObjects.twomers = new Array()
  pValueObjects.threemers = new Array()
  pValueObjects.fourmers = new Array()
  var i

  //parse onemers
  for (i = 0; i < parsed1.onemers.length; i++) {
    var onemer1 = parsed1.onemers[i]
    var n1 = 0
    Object.keys(onemer1).forEach(function(nucleotide) {
      n1 = n1 + onemer1[nucleotide].count
    })

    var onemer2 = parsed2.onemers[i]
    var n2 = 0
    Object.keys(onemer2).forEach(function(nucleotide) {
      n2 = n2 + onemer2[nucleotide].count
    })

    nucleotides = new Array()

    Object.keys(onemer1).forEach(function(nucleotide) {
      var p1 = onemer1[nucleotide].proportion
      var p2 = onemer2[nucleotide].proportion
      var y1 = onemer1[nucleotide].count
      var y2 = onemer2[nucleotide].count
      var po = (y1 + y2) / (n1 + n2)
      var z = (p1 - p2 - 0) / Math.sqrt(po * (1 - po) * (1 / n1 + 1 / n2))
      var pValue = GetZPercent(z)
      nucleotides.push(pValue)
      // console.log("ZSCORE: ", z)
    })

    var newOnemer = new Object()
    newOnemer = {
      A: nucleotides[0],
      C: nucleotides[1],
      G: nucleotides[2],
      T: nucleotides[3]
    }
    pValueObjects.onemers.push(newOnemer)
  }

  //parse twomers
  for (i = 0; i < parsed1.twomers.length; i++) {
    var twomer1 = parsed1.twomers[i]
    var n1 = 0
    Object.keys(twomer1).forEach(function(nucleotide) {
      n1 = n1 + twomer1[nucleotide].count
    })

    var twomer2 = parsed2.twomers[i]
    var n2 = 0
    Object.keys(twomer2).forEach(function(nucleotide) {
      n2 = n2 + twomer2[nucleotide].count
    })

    nucleotides = new Array()

    Object.keys(twomer1).forEach(function(nucleotide) {
      var p1 = twomer1[nucleotide].proportion
      var p2 = twomer2[nucleotide].proportion
      var y1 = twomer1[nucleotide].count
      var y2 = twomer2[nucleotide].count
      var po = (y1 + y2) / (n1 + n2)
      var z = (p1 - p2 - 0) / Math.sqrt(po * (1 - po) * (1 / n1 + 1 / n2))
      var pValue = GetZPercent(z)
      nucleotides.push(pValue)
    })
    var newTwomer = new Object()
    newTowmer = {
      A: nucleotides[0],
      C: nucleotides[1],
      G: nucleotides[2],
      T: nucleotides[3]
    }
    pValueObjects.twomers.push(newTwomer)
  }
  //parse threemers
  for (i = 0; i < parsed1.threemers.length; i++) {
    var threemer1 = parsed1.threemers[i]
    var n1 = 0
    Object.keys(threemer1).forEach(function(nucleotide) {
      n1 = n1 + threemer1[nucleotide].count
    })

    var threemer2 = parsed2.threemers[i]
    var n2 = 0
    Object.keys(threemer2).forEach(function(nucleotide) {
      n2 = n2 + threemer2[nucleotide].count
    })

    nucleotides = new Array()

    Object.keys(threemer1).forEach(function(nucleotide) {
      var p1 = threemer1[nucleotide].proportion
      var p2 = threemer2[nucleotide].proportion
      var y1 = threemer1[nucleotide].count
      var y2 = threemer2[nucleotide].count
      var po = (y1 + y2) / (n1 + n2)
      var z = (p1 - p2 - 0) / Math.sqrt(po * (1 - po) * (1 / n1 + 1 / n2))
      var pValue = GetZPercent(z)
      nucleotides.push(pValue)
    })

    var newThreemer = new Object()
    newThreemer = {
      A: nucleotides[0],
      C: nucleotides[1],
      G: nucleotides[2],
      T: nucleotides[3]
    }
    pValueObjects.threemers.push(newThreemer)
  }

  //parse fourmers
  for (i = 0; i < parsed1.fourmers.length; i++) {
    var fourmer1 = parsed1.fourmers[i]
    var n1 = 0
    Object.keys(fourmer1).forEach(function(nucleotide) {
      n1 = n1 + fourmer1[nucleotide].count
    })

    var fourmer2 = parsed2.fourmers[i]
    var n2 = 0
    Object.keys(fourmer2).forEach(function(nucleotide) {
      n2 = n2 + fourmer2[nucleotide].count
    })

    nucleotides = new Array()

    Object.keys(fourmer1).forEach(function(nucleotide) {
      var p1 = fourmer1[nucleotide].proportion
      var p2 = fourmer2[nucleotide].proportion
      var y1 = fourmer1[nucleotide].count
      var y2 = fourmer2[nucleotide].count
      var po = (y1 + y2) / (n1 + n2)
      var z = (p1 - p2 - 0) / Math.sqrt(po * (1 - po) * (1 / n1 + 1 / n2))
      var pValue = GetZPercent(z)
      nucleotides.push(pValue)
    })

    var newFourmer = new Object()
    newFourmer = {
      A: nucleotides[0],
      C: nucleotides[1],
      G: nucleotides[2],
      T: nucleotides[3]
    }
    pValueObjects.fourmers.push(newFourmer)
  }
  return pValueObjects
}

const GetZPercent = z => {
  //z == number of standard deviations from the mean

  //if z is greater than 6.5 standard deviations from the mean
  //the number of significant digits will be outside of a reasonable
  //range
  if (z < -6.5) return 0.0
  if (z > 6.5)
    //return 1.0;
    return 0.0

  var factK = 1
  var sum = 0
  var term = 1
  var k = 0
  var loopStop = Math.exp(-23)
  while (Math.abs(term) > loopStop) {
    term =
      0.3989422804 *
      Math.pow(-1, k) *
      Math.pow(z, k) /
      (2 * k + 1) /
      Math.pow(2, k) *
      Math.pow(z, k + 1) /
      factK
    sum += term
    k++
    factK *= k
  }
  sum += 0.5

  if (z < 0) {
    return sum * 2
  } else {
    return (1 - sum) * 2
  }
}

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
  console.log(pValueTable)
  for (let key in pValueTable) {
    let tableString = `<div>${key}</div>`
    tableString += '<table style="margin: 16px">'
    for (let j = 0; j < pValueTable[key].length; j++) {
      tableString += "<tr>"
      for (let nuc in pValueTable[key][j]) {
        tableString += `<td>${pValueTable[key][j][nuc]}</td>`
      }
      tableString += "</tr>"
    }
    tableString += "</table>"
    const renderNode = document.getElementById(`table${key}`)
    renderNode.innerHTML = tableString
  }
}

// Init function
;(function() {
  const goButton = document.getElementById("go")
  goButton.addEventListener("click", beginProcessing)
})()
