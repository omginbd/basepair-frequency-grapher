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
    newTwomer = {
      AA: nucleotides[0],
      AC: nucleotides[1],
      AG: nucleotides[2],
      AT: nucleotides[3],
      CA: nucleotides[4],
      CC: nucleotides[5],
      CG: nucleotides[6],
      CT: nucleotides[7], 
      GA: nucleotides[8],
      GC: nucleotides[9],
      GG: nucleotides[10],
      GT: nucleotides[11],
      TA: nucleotides[12],
      TC: nucleotides[13],
      TG: nucleotides[14],
      TT: nucleotides[15]
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
    position = 0
    newThreemer = {
      AAA: nucleotides[position++],
      AAC: nucleotides[position++],
      AAG: nucleotides[position++],
      AAT: nucleotides[position++],
      ACA: nucleotides[position++],
      ACC: nucleotides[position++],
      ACG: nucleotides[position++],
      ACT: nucleotides[position++],
      AGA: nucleotides[position++],
      AGC: nucleotides[position++],
      AGG: nucleotides[position++],
      AGT: nucleotides[position++],
      ATA: nucleotides[position++],
      ATC: nucleotides[position++],
      ATG: nucleotides[position++],
      ATT: nucleotides[position++],
      CAA: nucleotides[position++],
      CAC: nucleotides[position++],
      CAG: nucleotides[position++],
      CAT: nucleotides[position++],
      CCA: nucleotides[position++],
      CCC: nucleotides[position++],
      CCG: nucleotides[position++],
      CCT: nucleotides[position++],
      CGA: nucleotides[position++],
      CGC: nucleotides[position++],
      CGG: nucleotides[position++],
      CGT: nucleotides[position++],
      CTA: nucleotides[position++],
      CTC: nucleotides[position++],
      CTG: nucleotides[position++],
      CTT: nucleotides[position++],
      GAA: nucleotides[position++],
      GAC: nucleotides[position++],
      GAG: nucleotides[position++],
      GAT: nucleotides[position++],
      GCA: nucleotides[position++],
      GCC: nucleotides[position++],
      GCG: nucleotides[position++],
      GCT: nucleotides[position++],
      GGA: nucleotides[position++],
      GGC: nucleotides[position++],
      GGG: nucleotides[position++],
      GGT: nucleotides[position++],
      GTA: nucleotides[position++],
      GTC: nucleotides[position++],
      GTG: nucleotides[position++],
      GTT: nucleotides[position++],
      TAA: nucleotides[position++],
      TAC: nucleotides[position++],
      TAG: nucleotides[position++],
      TAT: nucleotides[position++],
      TCA: nucleotides[position++],
      TCC: nucleotides[position++],
      TCG: nucleotides[position++],
      TCT: nucleotides[position++],
      TGA: nucleotides[position++],
      TGC: nucleotides[position++],
      TGG: nucleotides[position++],
      TGT: nucleotides[position++],
      TTA: nucleotides[position++],
      TTC: nucleotides[position++],
      TTG: nucleotides[position++],
      TTT: nucleotides[position++]
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

    // Object.keys(fourmer1).forEach(function(nucleotide) {
    //   var p1 = fourmer1[nucleotide].proportion
    //   var p2 = fourmer2[nucleotide].proportion
    //   var y1 = fourmer1[nucleotide].count
    //   var y2 = fourmer2[nucleotide].count
    //   var po = (y1 + y2) / (n1 + n2)
    //   var z = (p1 - p2 - 0) / Math.sqrt(po * (1 - po) * (1 / n1 + 1 / n2))
    //   var pValue = GetZPercent(z)
    //   nucleotides.push(pValue)
    // })

    for (i = 0; i < 10; i++){
      nucleotides.push(i)
    }

    var newFourmer = new Object()
    var position = 0

    newFourmer = {
      AAAA: nucleotides[position++],
      AAAC: nucleotides[position++],
      AAAG: nucleotides[position++],
      AAAT: nucleotides[position++],
      AACA: nucleotides[position++],
      AACC: nucleotides[position++],
      AACG: nucleotides[position++],
      AACT: nucleotides[position++],
      AAGA: nucleotides[position++],
      AAGC: nucleotides[position++],
      AAGG: nucleotides[position++],
      AAGT: nucleotides[position++],
      AATA: nucleotides[position++],
      AATC: nucleotides[position++],
      AATG: nucleotides[position++],
      AATT: nucleotides[position++],
      ACAA: nucleotides[position++],
      ACAC: nucleotides[position++],
      ACAG: nucleotides[position++],
      ACAT: nucleotides[position++],
      ACCA: nucleotides[position++],
      ACCC: nucleotides[position++],
      ACCG: nucleotides[position++],
      ACCT: nucleotides[position++],
      ACGA: nucleotides[position++],
      ACGC: nucleotides[position++],
      ACGG: nucleotides[position++],
      ACGT: nucleotides[position++],
      ACTA: nucleotides[position++],
      ACTC: nucleotides[position++],
      ACTG: nucleotides[position++],
      ACTT: nucleotides[position++],
      AGAA: nucleotides[position++],
      AGAC: nucleotides[position++],
      AGAG: nucleotides[position++],
      AGAT: nucleotides[position++],
      AGCA: nucleotides[position++],
      AGCC: nucleotides[position++],
      AGCG: nucleotides[position++],
      AGCT: nucleotides[position++],
      AGGA: nucleotides[position++],
      AGGC: nucleotides[position++],
      AGGG: nucleotides[position++],
      AGGT: nucleotides[position++],
      AGTA: nucleotides[position++],
      AGTC: nucleotides[position++],
      AGTG: nucleotides[position++],
      AGTT: nucleotides[position++],
      ATAA: nucleotides[position++],
      ATAC: nucleotides[position++],
      ATAG: nucleotides[position++],
      ATAT: nucleotides[position++],
      ATCA: nucleotides[position++],
      ATCC: nucleotides[position++],
      ATCG: nucleotides[position++],
      ATCT: nucleotides[position++],
      ATGA: nucleotides[position++],
      ATGC: nucleotides[position++],
      ATGG: nucleotides[position++],
      ATGT: nucleotides[position++],
      ATTA: nucleotides[position++],
      ATTC: nucleotides[position++],
      ATTG: nucleotides[position++],
      ATTT: nucleotides[position++],	
      CAAA: nucleotides[position++],
      CAAC: nucleotides[position++],
      CAAG: nucleotides[position++],
      CAAT: nucleotides[position++],
      CACA: nucleotides[position++],
      CACC: nucleotides[position++],
      CACG: nucleotides[position++],
      CACT: nucleotides[position++],
      CAGA: nucleotides[position++],
      CAGC: nucleotides[position++],
      CAGG: nucleotides[position++],
      CAGT: nucleotides[position++],
      CATA: nucleotides[position++],
      CATC: nucleotides[position++],
      CATG: nucleotides[position++],
      CATT: nucleotides[position++],
      CCAA: nucleotides[position++],
      CCAC: nucleotides[position++],
      CCAG: nucleotides[position++],
      CCAT: nucleotides[position++],
      CCCA: nucleotides[position++],
      CCCC: nucleotides[position++],
      CCCG: nucleotides[position++],
      CCCT: nucleotides[position++],
      CCGA: nucleotides[position++],
      CCGC: nucleotides[position++],
      CCGG: nucleotides[position++],
      CCGT: nucleotides[position++],
      CCTA: nucleotides[position++],
      CCTC: nucleotides[position++],
      CCTG: nucleotides[position++],
      CCTT: nucleotides[position++],
      CGAA: nucleotides[position++],
      CGAC: nucleotides[position++],
      CGAG: nucleotides[position++],
      CGAT: nucleotides[position++],
      CGCA: nucleotides[position++],
      CGCC: nucleotides[position++],
      CGCG: nucleotides[position++],
      CGCT: nucleotides[position++],
      CGGA: nucleotides[position++],
      CGGC: nucleotides[position++],
      CGGG: nucleotides[position++],
      CGGT: nucleotides[position++],
      CGTA: nucleotides[position++],
      CGTC: nucleotides[position++],
      CGTG: nucleotides[position++],
      CGTT: nucleotides[position++],
      CTAA: nucleotides[position++],
      CTAC: nucleotides[position++],
      CTAG: nucleotides[position++],
      CTAT: nucleotides[position++],
      CTCA: nucleotides[position++],
      CTCC: nucleotides[position++],
      CTCG: nucleotides[position++],
      CTCT: nucleotides[position++],
      CTGA: nucleotides[position++],
      CTGC: nucleotides[position++],
      CTGG: nucleotides[position++],
      CTGT: nucleotides[position++],
      CTTA: nucleotides[position++],
      CTTC: nucleotides[position++],
      CTTG: nucleotides[position++],
      CTTT: nucleotides[position++],	
      GAAA: nucleotides[position++],
      GAAC: nucleotides[position++],
      GAAG: nucleotides[position++],
      GAAT: nucleotides[position++],
      GACA: nucleotides[position++],
      GACC: nucleotides[position++],
      GACG: nucleotides[position++],
      GACT: nucleotides[position++],
      GAGA: nucleotides[position++],
      GAGC: nucleotides[position++],
      GAGG: nucleotides[position++],
      GAGT: nucleotides[position++],
      GATA: nucleotides[position++],
      GATC: nucleotides[position++],
      GATG: nucleotides[position++],
      GATT: nucleotides[position++],
      GCAA: nucleotides[position++],
      GCAC: nucleotides[position++],
      GCAG: nucleotides[position++],
      GCAT: nucleotides[position++],
      GCCA: nucleotides[position++],
      GCCC: nucleotides[position++],
      GCCG: nucleotides[position++],
      GCCT: nucleotides[position++],
      GCGA: nucleotides[position++],
      GCGC: nucleotides[position++],
      GCGG: nucleotides[position++],
      GCGT: nucleotides[position++],
      GCTA: nucleotides[position++],
      GCTC: nucleotides[position++],
      GCTG: nucleotides[position++],
      GCTT: nucleotides[position++],
      GGAA: nucleotides[position++],
      GGAC: nucleotides[position++],
      GGAG: nucleotides[position++],
      GGAT: nucleotides[position++],
      GGCA: nucleotides[position++],
      GGCC: nucleotides[position++],
      GGCG: nucleotides[position++],
      GGCT: nucleotides[position++],
      GGGA: nucleotides[position++],
      GGGC: nucleotides[position++],
      GGGG: nucleotides[position++],
      GGGT: nucleotides[position++],
      GGTA: nucleotides[position++],
      GGTC: nucleotides[position++],
      GGTG: nucleotides[position++],
      GGTT: nucleotides[position++],
      GTAA: nucleotides[position++],
      GTAC: nucleotides[position++],
      GTAG: nucleotides[position++],
      GTAT: nucleotides[position++],
      GTCA: nucleotides[position++],
      GTCC: nucleotides[position++],
      GTCG: nucleotides[position++],
      GTCT: nucleotides[position++],
      GTGA: nucleotides[position++],
      GTGC: nucleotides[position++],
      GTGG: nucleotides[position++],
      GTGT: nucleotides[position++],
      GTTA: nucleotides[position++],
      GTTC: nucleotides[position++],
      GTTG: nucleotides[position++],
      GTTT: nucleotides[position++],	
      TAAA: nucleotides[position++],
      TAAC: nucleotides[position++],
      TAAG: nucleotides[position++],
      TAAT: nucleotides[position++],
      TACA: nucleotides[position++],
      TACC: nucleotides[position++],
      TACG: nucleotides[position++],
      TACT: nucleotides[position++],
      TAGA: nucleotides[position++],
      TAGC: nucleotides[position++],
      TAGG: nucleotides[position++],
      TAGT: nucleotides[position++],
      TATA: nucleotides[position++],
      TATC: nucleotides[position++],
      TATG: nucleotides[position++],
      TATT: nucleotides[position++],
      TCAA: nucleotides[position++],
      TCAC: nucleotides[position++],
      TCAG: nucleotides[position++],
      TCAT: nucleotides[position++],
      TCCA: nucleotides[position++],
      TCCC: nucleotides[position++],
      TCCG: nucleotides[position++],
      TCCT: nucleotides[position++],
      TCGA: nucleotides[position++],
      TCGC: nucleotides[position++],
      TCGG: nucleotides[position++],
      TCGT: nucleotides[position++],
      TCTA: nucleotides[position++],
      TCTC: nucleotides[position++],
      TCTG: nucleotides[position++],
      TCTT: nucleotides[position++],
      TGAA: nucleotides[position++],
      TGAC: nucleotides[position++],
      TGAG: nucleotides[position++],
      TGAT: nucleotides[position++],
      TGCA: nucleotides[position++],
      TGCC: nucleotides[position++],
      TGCG: nucleotides[position++],
      TGCT: nucleotides[position++],
      TGGA: nucleotides[position++],
      TGGC: nucleotides[position++],
      TGGG: nucleotides[position++],
      TGGT: nucleotides[position++],
      TGTA: nucleotides[position++],
      TGTC: nucleotides[position++],
      TGTG: nucleotides[position++],
      TGTT: nucleotides[position++],
      TTAA: nucleotides[position++],
      TTAC: nucleotides[position++],
      TTAG: nucleotides[position++],
      TTAT: nucleotides[position++],
      TTCA: nucleotides[position++],
      TTCC: nucleotides[position++],
      TTCG: nucleotides[position++],
      TTCT: nucleotides[position++],
      TTGA: nucleotides[position++],
      TTGC: nucleotides[position++],
      TTGG: nucleotides[position++],
      TTGT: nucleotides[position++],
      TTTA: nucleotides[position++],
      TTTC: nucleotides[position++],
      TTTG: nucleotides[position++],
      TTTT: nucleotides[position++]
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
