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

// TODO
const normalDistributionZTest = (kmer1, kmer2, n1, n2) => {}

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
                n2
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

const getPValue = function(zscore) {
  zscore = parseFloat(zscore)

  if (isNaN(zscore)) {
    throw new TypeError("zscore is not a valid number")
  }

  var yZscore = -3.4
  var xZscore = 0.09

  if (zscore === 0) {
    return 0.5
  }

  if (zscore > 0) {
    if (zscore > 3.49) {
      return 1
    }

    zscore = Math.floor(zscore * 100) / 100
    yZscore = Math.floor(zscore * 10) / 10
    yZscore = -yZscore
  } else {
    if (zscore < -3.49) {
      return 0
    }

    zscore = Math.ceil(zscore * 100) / 100
    yZscore = Math.ceil(zscore * 10) / 10
  }
  xZscore = Math.abs(Math.round((zscore % yZscore) * 10000) / 10000)

  var z100 = isNaN(xZscore) ? Math.abs(zscore) : xZscore
  var z10 = yZscore === 0 ? "0.0" : yZscore.toFixed(1)
  var col = ZTABLE.z.indexOf(z100)
  var perc = ZTABLE[z10][col]

  if (zscore > 0) {
    perc = Math.round((1 - perc) * 10000) / 10000
  }

  return perc
}

const ZTABLE = {
  z: [0.09, 0.08, 0.07, 0.06, 0.05, 0.04, 0.03, 0.02, 0.01, 0],
  "-3.4": [
    0.0002,
    0.0003,
    0.0003,
    0.0003,
    0.0003,
    0.0003,
    0.0003,
    0.0003,
    0.0003,
    0.0003
  ],
  "-3.3": [
    0.0003,
    0.0004,
    0.0004,
    0.0004,
    0.0004,
    0.0004,
    0.0004,
    0.0005,
    0.0005,
    0.0005
  ],
  "-3.2": [
    0.0005,
    0.0005,
    0.0005,
    0.0006,
    0.0006,
    0.0006,
    0.0006,
    0.0006,
    0.0007,
    0.0007
  ],
  "-3.1": [
    0.0007,
    0.0007,
    0.0008,
    0.0008,
    0.0008,
    0.0008,
    0.0009,
    0.0009,
    0.0009,
    0.001
  ],
  "-3.0": [
    0.001,
    0.001,
    0.0011,
    0.0011,
    0.0011,
    0.0012,
    0.0012,
    0.0013,
    0.0013,
    0.0013
  ],
  "-2.9": [
    0.0014,
    0.0014,
    0.0015,
    0.0015,
    0.0016,
    0.0016,
    0.0017,
    0.0018,
    0.0018,
    0.0019
  ],
  "-2.8": [
    0.0019,
    0.002,
    0.0021,
    0.0021,
    0.0022,
    0.0023,
    0.0023,
    0.0024,
    0.0025,
    0.0026
  ],
  "-2.7": [
    0.0026,
    0.0027,
    0.0028,
    0.0029,
    0.003,
    0.0031,
    0.0032,
    0.0033,
    0.0034,
    0.0035
  ],
  "-2.6": [
    0.0036,
    0.0037,
    0.0038,
    0.0039,
    0.004,
    0.0041,
    0.0043,
    0.0044,
    0.0045,
    0.0047
  ],
  "-2.5": [
    0.0048,
    0.0049,
    0.0051,
    0.0052,
    0.0054,
    0.0055,
    0.0057,
    0.0059,
    0.006,
    0.0062
  ],
  "-2.4": [
    0.0064,
    0.0066,
    0.0068,
    0.0069,
    0.0071,
    0.0073,
    0.0075,
    0.0078,
    0.008,
    0.0082
  ],
  "-2.3": [
    0.0084,
    0.0087,
    0.0089,
    0.0091,
    0.0094,
    0.0096,
    0.0099,
    0.0102,
    0.0104,
    0.0107
  ],
  "-2.2": [
    0.011,
    0.0113,
    0.0116,
    0.0119,
    0.0122,
    0.0125,
    0.0129,
    0.0132,
    0.0136,
    0.0139
  ],
  "-2.1": [
    0.0143,
    0.0146,
    0.015,
    0.0154,
    0.0158,
    0.0162,
    0.0166,
    0.017,
    0.0174,
    0.0179
  ],
  "-2.0": [
    0.0183,
    0.0188,
    0.0192,
    0.0197,
    0.0202,
    0.0207,
    0.0212,
    0.0217,
    0.0222,
    0.0228
  ],
  "-1.9": [
    0.0233,
    0.0239,
    0.0244,
    0.025,
    0.0256,
    0.0262,
    0.0268,
    0.0274,
    0.0281,
    0.0287
  ],
  "-1.8": [
    0.0294,
    0.0301,
    0.0307,
    0.0314,
    0.0322,
    0.0329,
    0.0336,
    0.0344,
    0.0351,
    0.0359
  ],
  "-1.7": [
    0.0367,
    0.0375,
    0.0384,
    0.0392,
    0.0401,
    0.0409,
    0.0418,
    0.0427,
    0.0436,
    0.0446
  ],
  "-1.6": [
    0.0455,
    0.0465,
    0.0475,
    0.0485,
    0.0495,
    0.0505,
    0.0516,
    0.0526,
    0.0537,
    0.0548
  ],
  "-1.5": [
    0.0559,
    0.0571,
    0.0582,
    0.0594,
    0.0606,
    0.0618,
    0.063,
    0.0643,
    0.0655,
    0.0668
  ],
  "-1.4": [
    0.0681,
    0.0694,
    0.0708,
    0.0721,
    0.0735,
    0.0749,
    0.0764,
    0.0778,
    0.0793,
    0.0808
  ],
  "-1.3": [
    0.0823,
    0.0838,
    0.0853,
    0.0869,
    0.0885,
    0.0901,
    0.0918,
    0.0934,
    0.0951,
    0.0968
  ],
  "-1.2": [
    0.0985,
    0.1003,
    0.102,
    0.1038,
    0.1056,
    0.1075,
    0.1093,
    0.1112,
    0.1131,
    0.1151
  ],
  "-1.1": [
    0.117,
    0.119,
    0.121,
    0.123,
    0.1251,
    0.1271,
    0.1292,
    0.1314,
    0.1335,
    0.1357
  ],
  "-1.0": [
    0.1379,
    0.1401,
    0.1423,
    0.1446,
    0.1469,
    0.1492,
    0.1515,
    0.1539,
    0.1562,
    0.1587
  ],
  "-0.9": [
    0.1611,
    0.1635,
    0.166,
    0.1685,
    0.1711,
    0.1736,
    0.1762,
    0.1788,
    0.1814,
    0.1841
  ],
  "-0.8": [
    0.1867,
    0.1894,
    0.1922,
    0.1949,
    0.1977,
    0.2005,
    0.2033,
    0.2061,
    0.209,
    0.2119
  ],
  "-0.7": [
    0.2148,
    0.2177,
    0.2206,
    0.2236,
    0.2266,
    0.2296,
    0.2327,
    0.2358,
    0.2389,
    0.242
  ],
  "-0.6": [
    0.2451,
    0.2483,
    0.2514,
    0.2546,
    0.2578,
    0.2611,
    0.2643,
    0.2676,
    0.2709,
    0.2743
  ],
  "-0.5": [
    0.2776,
    0.281,
    0.2843,
    0.2877,
    0.2912,
    0.2946,
    0.2981,
    0.3015,
    0.305,
    0.3085
  ],
  "-0.4": [
    0.3121,
    0.3156,
    0.3192,
    0.3228,
    0.3264,
    0.33,
    0.3336,
    0.3372,
    0.3409,
    0.3446
  ],
  "-0.3": [
    0.3483,
    0.352,
    0.3557,
    0.3594,
    0.3632,
    0.3669,
    0.3707,
    0.3745,
    0.3783,
    0.3821
  ],
  "-0.2": [
    0.3829,
    0.3897,
    0.3936,
    0.3974,
    0.4013,
    0.4052,
    0.409,
    0.4129,
    0.4168,
    0.4207
  ],
  "-0.1": [
    0.4247,
    0.4286,
    0.4325,
    0.4364,
    0.4404,
    0.4443,
    0.4483,
    0.4522,
    0.4562,
    0.4602
  ],
  "0.0": [
    0.4641,
    0.4681,
    0.4721,
    0.4761,
    0.4801,
    0.484,
    0.488,
    0.492,
    0.496,
    0.5
  ]
}
