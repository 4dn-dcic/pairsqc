// This function requires d3.v3 and provides a simple template for plotting an interactive multi-line graph. 
// When mouse is over one of the lines, the line changes color and its label is displayed.
//
// Usage example:
//
// tsvcolfile = 'tsvcol_for_somedatafile.tsv';
// d3.tsv(tsvcolfile, function(tsvcolumns) {
//  tsvfile = 'somedatafile.tsv';
//  interactive_multiline_plot(tsvfile, tsvcolumns, 1, 9, 0, 4000, 'log10 distance', 'read count', 'my_div_id');
// });
//
// A tsvcolfile specifies individual lines as columnx and columny, matching the x, y column names in a separately given multi-column tsvfile.
// It also defines on- and off- colors (color and offcolor) and labels (name).
// An example tsvcolfile looks as below (tab-delimited):
//
// columnx columny color   offcolor        name
// distance        count.Inner     #ff0000 #999999 Inner
// distance        count.Outer     #00ff00 #999999 Outer
// ...
//
// Author: Soo Lee (duplexa@gmail.com), Department of Biomedical Informatics, Harvard Medical School.

function interactive_multiline_plot(tsvfile, tsvcolumns, xmin, xmax, ymin, ymax, x_label, y_label, div_id){
    
    // canvas size. (data-independent)
    var WIDTH = 300; 
    var HEIGHT = 300;
    var MARGINS = {
        top: 20,
        right: 20,
        bottom: 70,
        left: 100
    };
    var width = WIDTH + MARGINS.left + MARGINS.right;
    var height = HEIGHT + MARGINS.top + MARGINS.bottom; 
    
    // X and Y scale (could be data-dependent)
    var xScale = d3.scale.linear().range([MARGINS.left, width - MARGINS.right]).domain([xmin,xmax]);
    var yScale = d3.scale.linear().range([height - MARGINS.bottom, MARGINS.top]).domain([ymin,ymax]);
    
    // start (data-independent)
    var vis = d3.select('#' + div_id).append("svg")
                .attr("width", width)
                .attr("height", height);
    
    // Axes (data-independent)
    xAxis = d3.svg.axis()
        .scale(xScale),
    vis.append("svg:g")
        .attr("class","__axis__") 
        .attr("fill", "none")
        .attr("stroke","#000")
        .attr("stroke-width", 1)
        .attr("shape-rendering", "crispEdges")
        .attr("transform", "translate(0," + (HEIGHT + MARGINS.top) + ")")
        .call(xAxis);
    yAxis = d3.svg.axis()
        .scale(yScale)
        .orient("left");
    vis.append("svg:g")
        .attr("class","__axis__")  
        .attr("fill", "none")
        .attr("stroke","#000")
        .attr("stroke-width", 1)
        .attr("shape-rendering", "crispEdges")
        .attr("transform", "translate(" + (MARGINS.left) + ",0)")
        .call(yAxis);
    vis.selectAll(".__axis__ text")
        .style("font-family", "Lato")
        .style("font-size", "13px")
        .attr("fill", "#000");

    // Axis labels (data-independent)
    vis.append("text")
       .attr("transform", "translate(" + (WIDTH/2 + MARGINS.left) + " ," + (height - MARGINS.bottom/3) + ")")
       .style("text-anchor", "middle")
       .text(x_label);
    vis.append("text")
        .attr({y:MARGINS.left/3, x:-HEIGHT/2})
        .attr("transform", "rotate(-90)")
        .style("text-anchor", "middle")
        .text(y_label);

    // data to line transformation function (data-independent)
    var lineGen = d3.svg.line()
      .x(function(d) {
        return xScale(d.x);
      })
      .y(function(d) {
        return yScale(d.y);
      })
      .interpolate("basis");

    // clip-area for plot lines, so that they don't go beyond plot region.
    vis.append("clipPath")
       .attr("id", div_id + "__clip_rect__")
       .append("svg:rect")
       .attr({x:MARGINS.left, y:MARGINS.top})
       .attr("width", WIDTH)
       .attr("height", HEIGHT);

    // plotting data in various colors (data-independent)
    function plot_default(){
      tsvcolumns.forEach(function(t){
      
        var field_key_y = t.columny;
        var field_key_x = t.columnx;
        var lineColor = t.color;  // color when the column is on
        var offColor = t.offcolor;  // default color when the column is off
        var lineName = t.name; // label
      
        d3.tsv(tsvfile, function(data) {
          data.forEach(function(d) {
             d.x = d[field_key_x] = +d[field_key_x];
             d.y = d[field_key_y] = +d[field_key_y];
          });
      
          // adding data lines (data-dependent)
          vis.append('svg:path')
            .attr('id', '__path__' + field_key_y)
            .attr('d', lineGen(data))
            .attr('stroke', offColor)  // start with default color
            .attr('stroke-width', 1)
            .attr('fill', 'none')
            .attr('clip-path', 'url(#' + div_id + '__clip_rect__)')
            .on("mouseover", function(){
               var p1 = d3.mouse(this);
               d3.select(this).attr('stroke', lineColor).attr('stroke-width', 3);
               var label_class = '__label__';
               vis.append("text").attr('class',label_class).attr({x:p1[0]+10, y:p1[1]-10}).text(lineName).style("font-size", "20px");
            })
            .on("mouseout", function(){
               d3.select(this).attr('stroke', offColor).attr('stroke-width', 1);
               var label_class = '.__label__';
               d3.selectAll(label_class).remove();
            });
     
        }); 
      });
    }
    
    plot_default();
}


