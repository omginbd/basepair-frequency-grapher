const tables = { onemers: "", twomers: "", threemers: "", fourmers: "" }

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
  } else if (z == 0){
    return 0
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
    const sigLevel = document.querySelector("#sigLevel").value
    const testTypes = document.getElementsByName("test-type")
    let test
    for (let type of testTypes) test = type.checked ? type.value : test
    updateTables(
      compareReads(parseFile(f1.result), parseFile(f2.result), test),
      sigLevel
    )
  } else {
    setTimeout(awaitProcessing, 100, f1, f2)
  }
}

const twoProportionZTest = (kmer1, kmer2, n1, n2) => {
  const { count: y1, proportion: p1 } = kmer1
  const { count: y2, proportion: p2 } = kmer2
  const po = (y1 + y2) / (n1 + n2)
  const z = (p1 - p2) / Math.sqrt(po * (1 - po) * (1 / n1 + 1 / n2))
  return GetZPercent(z)
}

// Source: http://www.math.ucla.edu/~tom/distributions/normal.html

const normalcdf = (X) => {   //HASTINGS.  MAX ERROR = .000001
	var T=1/(1+.2316419*Math.abs(X))
	var D=.3989423*Math.exp(-X*X/2)
	var Prob=D*T*(.3193815+T*(-.3565638+T*(1.781478+T*(-1.821256+T*1.330274))))
	if (X>0) {
		Prob=1-Prob
	}
	return Prob
}

const compute = (Z, M, SD) => {
  with (Math) {
    if (SD<0) {
      alert("The standard deviation must be nonnegative.")
    } else if (SD==0) {
        if (Z<M){
            Prob=0
        } else {
          Prob=1
      }
    } else {
      Prob=normalcdf((Z-M)/SD)
      Prob=round(100000*Prob)/100000
    }
  }
  return Prob
}

const normalDistributionZTest = (kmer1, kmer2, n1, n2, divisor) => {
  const { count: y1, proportion: p1 } = kmer1
  const { count: y2, proportion: p2 } = kmer2
  const sd = Math.sqrt(divisor * (p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2))
  const m = p1 - p2
  const znorm = m / sd

  const v = compute(znorm, m, sd)

  if (znorm < 0) {
    return v * 2
  } else if (znorm > 0) {
    return (1 - v) * 2
  } else {
    return 0
  }
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
const compareReads = (file1, file2, test) => {
  const divisor = document.querySelector("#divisor").value
  const toReturn = { onemers: [], twomers: [], threemers: [], fourmers: [] }
  for (let k in file1) {
    const collection = [file1[k], file2[k]]
    for (let [position] of collection[0].entries()) {
      const pValueObj = {}
      const n1 = Object.keys(collection[0][position]).reduce(
        (memo, val) => memo + collection[0][position][val].count,
        0
      )
      const n2 = Object.keys(collection[1][position]).reduce(
        (memo, val) => memo + collection[1][position][val].count,
        0
      )
      for (let kmer of Object.keys(collection[0][position])) {
        pValueObj[kmer] =
          test === "twoProportionZTest"
            ? twoProportionZTest(
                collection[0][position][kmer],
                collection[1][position][kmer],
                n1,
                n2
              )
            : normalDistributionZTest(
                collection[0][position][kmer],
                collection[1][position][kmer],
                n1,
                n2,
                divisor
              )
      }
      toReturn[k].push(pValueObj)
    }
  }
  return toReturn
}

const updateTables = (pValueTable, sigLevel) => {
  for (let key in pValueTable) {
    let tableString = `<div class="label">${key}</div>`
    tableString += '<table style="margin: 16px; display: block; width: 99vw">'
    tableString += "<tr><th />"
    for (let j = 0; j < Object.keys(pValueTable[key][0]).length; j++) {
      tableString += `<th class="theader">${
        Object.keys(pValueTable[key][0])[j]
      }</th>`
    }
    tableString += "</tr>"
    for (let j = 0; j < pValueTable[key].length; j++) {
      tableString += `<tr><td>${j}</td>`
      for (let nuc in pValueTable[key][j]) {
        const val = pValueTable[key][j][nuc]
        tableString += `<td class=${
          val < sigLevel ? "" : "insignificant"
        }>${val.toFixed(5)}</td>`
      }
      tableString += "</tr>"
    }
    tableString += "</table>"
    tables[key] = tableString
  }
  document.querySelector("#buttonWrapper").style.display = "flex"
  document.querySelector("#go").style.display = "none"
}

const changeTable = key => {
  const renderNode = document.querySelector("#table")
  renderNode.innerHTML = tables[key]
}

// Init function
;(function() {
  const goButton = document.getElementById("go")
  const k1Button = document.getElementById("onemers-button")
  const k2Button = document.getElementById("twomers-button")
  const k3Button = document.getElementById("threemers-button")
  const k4Button = document.getElementById("fourmers-button")
  goButton.addEventListener("click", beginProcessing)
  k1Button.addEventListener("click", () => changeTable("onemers"))
  k2Button.addEventListener("click", () => changeTable("twomers"))
  k3Button.addEventListener("click", () => changeTable("threemers"))
  k4Button.addEventListener("click", () => changeTable("fourmers"))
})()
