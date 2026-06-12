#show figure: set block(breakable: true)
#figure( // start preamble figure
  caption: text([Empirical rejection rates over 1000 Monte Carlo replications at $alpha = .05$, by model-comparison method and number of bootstrap resamples $B$. The $delta = 0$ columns give the Type I error rate; $delta > 0$ columns give power.]),
  kind: table,
  supplement: "Table", // end preamble figure

block[ // start block

  #let style-dict = (
    // tinytable style-dict after
    "0_0": 0, "0_1": 0, "1_1": 0, "2_1": 0, "3_1": 0, "4_1": 0, "5_1": 0, "6_1": 0, "7_1": 0, "8_1": 0, "9_1": 0, "10_1": 0, "11_1": 0, "12_1": 0, "13_1": 0, "14_1": 0, "15_1": 0, "16_1": 0, "0_2": 0, "1_2": 0, "2_2": 0, "3_2": 0, "4_2": 0, "5_2": 0, "6_2": 0, "7_2": 0, "8_2": 0, "9_2": 0, "10_2": 0, "11_2": 0, "12_2": 0, "13_2": 0, "14_2": 0, "15_2": 0, "16_2": 0, "0_3": 0, "1_3": 0, "2_3": 0, "3_3": 0, "4_3": 0, "5_3": 0, "6_3": 0, "7_3": 0, "8_3": 0, "9_3": 0, "10_3": 0, "11_3": 0, "12_3": 0, "13_3": 0, "14_3": 0, "15_3": 0, "16_3": 0, "0_4": 0, "1_4": 0, "2_4": 0, "3_4": 0, "4_4": 0, "5_4": 0, "6_4": 0, "7_4": 0, "8_4": 0, "9_4": 0, "10_4": 0, "11_4": 0, "12_4": 0, "13_4": 0, "14_4": 0, "15_4": 0, "16_4": 0, "0_5": 0, "1_5": 0, "2_5": 0, "3_5": 0, "4_5": 0, "5_5": 0, "6_5": 0, "7_5": 0, "8_5": 0, "9_5": 0, "10_5": 0, "11_5": 0, "12_5": 0, "13_5": 0, "14_5": 0, "15_5": 0, "16_5": 0, "0_6": 0, "1_6": 0, "2_6": 0, "3_6": 0, "4_6": 0, "5_6": 0, "6_6": 0, "7_6": 0, "8_6": 0, "9_6": 0, "10_6": 0, "11_6": 0, "12_6": 0, "13_6": 0, "14_6": 0, "15_6": 0, "16_6": 0, "0_7": 0, "1_7": 0, "2_7": 0, "3_7": 0, "4_7": 0, "5_7": 0, "6_7": 0, "7_7": 0, "8_7": 0, "9_7": 0, "10_7": 0, "11_7": 0, "12_7": 0, "13_7": 0, "14_7": 0, "15_7": 0, "16_7": 0, "0_8": 0, "1_8": 0, "2_8": 0, "3_8": 0, "4_8": 0, "5_8": 0, "6_8": 0, "7_8": 0, "8_8": 0, "9_8": 0, "10_8": 0, "11_8": 0, "12_8": 0, "13_8": 0, "14_8": 0, "15_8": 0, "16_8": 0, "0_9": 0, "1_9": 0, "2_9": 0, "3_9": 0, "4_9": 0, "5_9": 0, "6_9": 0, "7_9": 0, "8_9": 0, "9_9": 0, "10_9": 0, "11_9": 0, "12_9": 0, "13_9": 0, "14_9": 0, "15_9": 0, "16_9": 0
  )

  #let style-array = ( 
    // tinytable cell style after
    (align: center,),
  )

  // Helper function to get cell style
  #let get-style(x, y) = {
    let key = str(y) + "_" + str(x)
    if key in style-dict { style-array.at(style-dict.at(key)) } else { none }
  }

  // tinytable align-default-array before
  #let align-default-array = ( left, left, left, left, left, left, left, left, left, left, ) // tinytable align-default-array here
  #show table.cell: it => {
    if style-array.len() == 0 { return it }
    
    let style = get-style(it.x, it.y)
    if style == none { return it }
    
    let tmp = it
    if ("fontsize" in style) { tmp = text(size: style.fontsize, tmp) }
    if ("color" in style) { tmp = text(fill: style.color, tmp) }
    if ("indent" in style) { tmp = pad(left: style.indent, tmp) }
    if ("underline" in style) { tmp = underline(tmp) }
    if ("italic" in style) { tmp = emph(tmp) }
    if ("bold" in style) { tmp = strong(tmp) }
    if ("mono" in style) { tmp = math.mono(tmp) }
    if ("strikeout" in style) { tmp = strike(tmp) }
    if ("smallcaps" in style) { tmp = smallcaps(tmp) }
    tmp
  }

  #align(center, [

  #table( // tinytable table start
    column-gutter: 5pt,
    columns: (auto, auto, auto, auto, auto, auto, auto, auto, auto, auto),
    stroke: none,
    rows: auto,
    align: (x, y) => {
      let style = get-style(x, y)
      if style != none and "align" in style { style.align } else { left }
    },
    fill: (x, y) => {
      let style = get-style(x, y)
      if style != none and "background" in style { style.background }
    },
 table.hline(y: 1, start: 2, end: 6, stroke: 0.03em + black), table.hline(y: 1, start: 6, end: 10, stroke: 0.03em + black),
 table.hline(y: 2, start: 0, end: 10, stroke: 0.05em + black),
 table.hline(y: 17, start: 0, end: 10, stroke: 0.08em + black),
 table.hline(y: 0, start: 0, end: 10, stroke: 0.08em + black),
    // tinytable lines before

    // tinytable header start
    table.header(
      repeat: true,
[ ], [ ], table.cell(colspan: 4, align: center)[R²], table.cell(colspan: 4, align: center)[adj-R²],
[Method], [$B$], [0], [0.1], [0.2], [0.4], [0], [0.1], [0.2], [0.4],
    ),
    // tinytable header end

    // tinytable cell content after
table.cell(colspan: 10)[$n_1\/n_2 = 100\/100$],
[H&T], [100], [1.000], [1.000], [1.000], [1.000], [0.874], [0.924], [0.979], [1.000],
[H&T], [1000], [1.000], [1.000], [1.000], [1.000], [0.975], [0.985], [0.997], [1.000],
[Permutation], [100], [0.038], [0.101], [0.327], [0.908], [0.038], [0.101], [0.327], [0.908],
[Permutation], [1000], [0.044], [0.121], [0.347], [0.916], [0.044], [0.121], [0.347], [0.916],
table.cell(colspan: 10)[$n_1\/n_2 = 140\/60$],
[H&T], [100], [1.000], [1.000], [1.000], [1.000], [0.887], [0.918], [0.967], [0.999],
[H&T], [1000], [1.000], [1.000], [1.000], [1.000], [0.967], [0.974], [0.992], [1.000],
[Permutation], [100], [0.052], [0.106], [0.279], [0.848], [0.052], [0.106], [0.279], [0.848],
[Permutation], [1000], [0.056], [0.095], [0.270], [0.849], [0.056], [0.095], [0.270], [0.849],
table.cell(colspan: 10)[$n_1\/n_2 = 180\/20$],
[H&T], [100], [1.000], [1.000], [1.000], [1.000], [0.864], [0.879], [0.895], [0.978],
[H&T], [1000], [1.000], [1.000], [1.000], [1.000], [0.954], [0.960], [0.969], [0.991],
[Permutation], [100], [0.044], [0.071], [0.141], [0.447], [0.044], [0.071], [0.141], [0.447],
[Permutation], [1000], [0.042], [0.092], [0.122], [0.435], [0.042], [0.092], [0.122], [0.435],

    // tinytable footer after

    table.footer(
      repeat: false,
      // tinytable notes after
    table.cell(align: left, colspan: 10, text([_Note._ Columns within each metric block index the slope difference $delta$. Cells are the proportion of replications rejecting $H_0$ at $alpha = .05$; $delta = 0$ is the Type I error rate and $delta > 0$ is power. Row groups give the part-1/part-2 sample sizes; the deviating subgroup is part 2.])),
    ),
    

  ) // end table

  ]) // end align

] // end block
) // end figure
<tbl-results>

