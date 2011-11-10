<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
    <xsl:output encoding="UTF-8" method="html" omit-xml-declaration="yes" indent="yes"/>

    <xsl:variable name="width" select="1900"/>

    <xsl:template match="prsm">
        <html>
            <title>Scan #<xsl:value-of select="scan/scan-id"/> matched to protein #<xsl:value-of select="protein/protein-id"/></title>
            <script type="text/javascript" src="../../../../js/prsm.js"/>
            <body>
                <h3>Scan #<xsl:value-of select="scan/scan-id"/></h3>
                <div>Precursor Mass <xsl:value-of select="scan/precursor-mass"/>, <xsl:value-of select="count(scan/peaks/peak)"/> peaks. </div>
                <h3>Protein #<xsl:value-of select="protein/protein-id"/>
                    <xsl:text> </xsl:text>
                    <xsl:value-of select="protein/protein-name"/></h3>
                <script>
                    var sequence = '<xsl:value-of select="protein/protein-sequence"/>';
                </script>

                <div>
                Scale <input id="scale" value="100" size="5" onchange="repaintPrsm();"/>%
                    <a href="#" onclick="zoom(1);" id="zoomIn" alt="Zoom In">+Zoom In</a>
                    <xsl:text> </xsl:text>
                    <a href="#" onclick="zoom(-1);" id="zoomOut" alt="Zoom Out">-Zoom Out</a>
                    <xsl:text> </xsl:text>
                    <a href="#" onclick="prefixLen++; leftRightMove();" id="leftArrow" alt="Move protein to the left"><xsl:text disable-output-escaping="yes">&amp;larr;Left</xsl:text></a>
                    <xsl:text> </xsl:text>
                    <a href="#" onclick="prefixLen--; leftRightMove();" id="rightArrow" alt="Move protein to the right"><xsl:text disable-output-escaping="yes">Right&amp;rarr;</xsl:text></a>
                </div>

                <div id="prefix">...</div>


                <canvas id="prsm" width="6000px" height="1000"/>


                <textarea id="deltas" cols="50" rows="5"></textarea>

                <script>
                    var prsmCanvas = document.getElementById('prsm');
                    var leftArrow = document.getElementById('leftArrow');
                    var rightArrow = document.getElementById('rightArrow');
                    var ctx = prsmCanvas.getContext('2d');
                    var font = "10pt Arial";
                    ctx.font = font;
                    var peaks = [
                        <xsl:apply-templates select="scan/peaks/peak"/>
                    ];

                    var scaleControl =  document.getElementById('scale');

                    var prefixLen = 0;
                    initPrsm();

                    document.getElementById('deltas').value = deltas;
                    var prefix = document.getElementById('prefix');

                    function doKeyDown(e) {
                        switch (e.keyCode) {
                            case 107: zoom(+1); break;
                            case 109: zoom(-1); break;
                            case 39: if (sequence.length > prefixLen) prefixLen++; leftRightMove(); break;
                            case 37: if (prefixLen > 0) prefixLen--; leftRightMove(); break;
                        }
                    }

                    function leftRightMove(){
                        prefix.innerHTML = sequence.substr(0, prefixLen) + "...";
                        leftArrow.style.display = sequence.length > prefixLen ? "" : "none";
                        rightArrow.style.display = prefixLen > 0 ? "" : "none" ;
                        repaintPrsm();
                    }

                    function zoom(delta){
                        scaleControl.value = parseInt(scaleControl.value) + parseInt(delta);
                        if (1 > scaleControl.value) {
                            scaleControl.value = 1;
                        }
                        document.getElementById('zoomOut').style.display = scaleControl.value > 1 ? "" : "none";
                        repaintPrsm();
                    }

                    leftRightMove();

                    window.addEventListener('keydown', doKeyDown, true);
                </script>
            </body>
        </html>
    </xsl:template>


    <xsl:template match="peak">
        {mass: <xsl:value-of select="mass"/>,
        intencity: <xsl:value-of select="intencity"/>,
        charge: <xsl:value-of  select="charge"/>},
    </xsl:template>


</xsl:stylesheet>