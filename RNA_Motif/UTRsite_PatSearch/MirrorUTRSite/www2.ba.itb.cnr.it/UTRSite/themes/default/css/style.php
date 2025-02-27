
/******************************************************************************/
/*                                  LAYOUT CSS                                */
/******************************************************************************/

body {
    margin: 0;
    padding: 0;
    font: 12px Trebuchet MS, Verdana, Arial, Helvetica, sans-serif;
    color: #000000;
}

input, select, textarea {
    font-family: Verdana, Arial, Helvetica, sans-serif;
    font-size: 10px;
}

/******************************* LAYOUT : HEADER ******************************/

#sgl-header {
    background-color: #000080;
    height: 50px;
    margin: 0;
    padding: 0;
}
#sgl-header-left {
    float: left;
    margin: 0;
    padding: 5px 0 5px 15px;
    color: #ffffff;
}
#sgl-header-left img {
    vertical-align: middle;
}
#sgl-header-right {
    float: right;
    padding: 10px 10px 0 0;
    text-align: right;
    font-size: 10px;
    color: #ffff00;
}
#sgl-header-right a {
    color: #ffffff;
    cursor: pointer;
    text-decoration: none;
}
#sgl-header-right a:hover {
    color: #00FFFF;
    cursor: pointer;
    text-decoration: none;
    text-decoration: underline;
}
#sgl-header-right img {
    margin-left: 10px;
    margin-right: 5px;
    width: 20px;
    height: 20px;
    vertical-align: middle;
}
#sgl-header-right #headerLogAction {
    float: left;
    padding: 2px;
}
#sgl-header-right #headerUserInfo {
    float: left;
    padding-top: 5px;
}
#sgl-header-right #headerUserInfo .guest {
    font-weight: bold;
}

/***************************** LAYOUT : NAVIGATION ****************************/

#sgl-nav {
    font-size: 10px;
}

/***************************** LAYOUT : TABLES ********************************/

table {
    border: none;
    /* This is not a typo, we want first set a fallback for IE, then set the
     * real margin for real browsers ;) */
    margin: 0 auto;
}
td {
    padding: 2px 4px;
    font-size: 12px;
}
th {
    background-color:#AABBCC;
    text-align: left;
    color: #FFFFFF;
    font-size: 13px;
    line-height: 20px;
    padding: 4px 10px 4px 10px;
    border: 2px outset;
    white-space: nowrap;
}
#imRead {
    background-color: #c6c6c6;
}
#moduleOverview {
    width: 240px;
    height: 125px;
}
.wide {
    width: 90%;
}
.narrow {
    width: 60%;
}
.sgl-row-light {
    background-color: #efefef;
}
.sgl-row-dark {
    background-color: #dddddd;
}
.sgl-row-light td {
    font: 11px Verdana, Arial, Helvetica, sans-serif;
    line-height: 20px;
    border: 1px outset;
}
.sgl-row-dark td {
    font: 11px Verdana, Arial, Helvetica, sans-serif;
    line-height: 20px;
    border: 1px outset;
}

/****************************** LAYOUT : MAIN *********************************/

#sgl-main {
    top: 110px;
    margin-top: 10px;
}

/************************ LAYOUT : LEFT & RIGHT BLOCKS ************************/

#sgl-blocks-left, #sgl-blocks-right {
    position: absolute;
    margin-top: 110px;
    top: 0;
    z-index:1;
}
#sgl-blocks-left {
    width: 210px;
    left: 0;
}
#sgl-blocks-right {
    width: 180px;
    right: 0;
}
.options-block {
    margin: 20px 0;
}
div.sgl-blocks-left-item {
    padding: 4px 0 0 1px;
    margin: 0 5px 10px 5px;
}
div.sgl-blocks-right-item {
    padding: 4px 1px 0 0;
    margin: 0 5px 10px 5px;
}
.sgl-blocks-left-item-title, .sgl-blocks-right-item-title {
    background-color: #000080;
    background-image: url(bg_block.png);
    background-repeat: repeat-x;
    font-size: 14px;
    line-height: 22px;
    font-weight: bold;
    color: #FFFF80;
    text-align: center;
    border: 2px outset #efefef;
}
div.sgl-blocks-left-item-body, div.sgl-blocks-right-item-body {
    background-color: #EFF3F8;
    color: #0206AE;
    font-size: 11px;
    padding: 5px 10px;
    border-right: 2px outset #efefef;
    border-left: 2px outset #efefef;
    border-bottom: 2px outset #efefef;
}

/*************************** LAYOUT : MIDDLE BLOCKS ***************************/

#sgl-blocks-middle, #sgl-blocks-middle-nocols, #sgl-blocks-middle-leftcol, #sgl-blocks-middle-rightcol {
    position: relative;
    margin: 0 180px 0 210px;
    width: auto;
    min-width: 20%;
    font-size: 11px;
    z-index: 2;
    padding: 0 20px;
}
#sgl-blocks-middle #options {
    float: right;
    width: 28%;
}
#sgl-blocks-middle-nocols {
    margin: 0;
    text-align: center;
}
#sgl-blocks-middle-leftcol {
    margin: 0 0 0 210px;
}
#sgl-blocks-middle-rightcol {
    margin: 0 180px 0 0;
}
div.sgl-blocks-middle-title {
    background-image: url(bg_block2.png);
    background-repeat: repeat-x;
    background-color: #000080;
    color: #FFFF80;
    font-size: 14px;
    line-height: 22px;
    font-weight: bold;
    text-align: center;
    border: 2px outset #efefef;
}
div.sgl-blocks-middle-body {
    background-color: #ffffff;
    background-color: #EFF3F8;
    color: #0206AE;
    border: 1px solid #000080;
    border-top: none;
    text-align: center;
    border: 2px outset #efefef;
}

/* Holly Hack here so that tooltips don't act screwy:
 * http://www.positioniseverything.net/explorer/threepxtest.html */
/* Hide next from Mac IE plus non-IE \*/
* html #sgl-blocks-middle {
    height: 1%;
}
/* End hide from IE5/mac plus non-IE */

/******************************* LAYOUT : FOOTER ******************************/

#sgl-footer {
    position: relative;
    float: middle;
    clear: both;
    margin-bottom: 5px;
    padding-top: 10px;
    padding-bottom: 10px;
    font-size: 10px;
    text-align: center;
}

/******************************************************************************/
/*                                 CONTENT CSS                                */
/******************************************************************************/

/***************************** CONTENT : HEADINGS *****************************/
h1 {
    font-size: 26px;
    font-weight: bold;
    text-align: center;
}
h2 {
    font-size: 18px;
    text-align: center;
}
h3 {
    font-size: 16px;
    text-align: center;
}
h4 {
    font-size: 14px;
}
.Title {
    color: #0206AE;
    font-size: 22px;
}
.pageTitle {
    color: #0206AE;
    font-size: 22px;
    font-weight: normal;
}
.homeTitle {
    font-weight: bold;
    font-size: 14px;
}
.homeBody {
    font-size: 12px;
}

/***************************** CONTENT : ANCHORS ******************************/

a, a:visited {
    color: #0206AE;
    font-weight: bold;
    text-decoration: none;
}
a:hover {
    color: #9dcdfe;
}

/******************************* CONTENT : BLOCKS *****************************/

img.blocksAvatar {
    /* move the image up to be flush with bottom of title */
    position: relative;
    top: -5px;
    float: right;
    padding-left: 5px;
    align: left;
}

/*************************** CONTENT : MISCELLANEOUS **************************/

acronym {
    cursor: help;
}
hr {
    border: none;
    border: 2px dashed #999999;
}
img {
    border: none;
}
.alignCenter {
    text-align: center;
}
.error {
    color: #ff0000;
    font-size:12px;
    background: transparent;
}
.small {
    font-size: 0.8em;
}
.strong {
    font-weight: bold;
}
.title {
    color: #999999;
    font-weight: normal;
    font-size: 1.5em;
}
.detail {
    color: #999999;
    font-weight: normal;
    font-size: 0.8em;
}
.toolBtnSeparate {
    margin-left: 20px;
}
/*************************** MODULE: PUBLISHER ********************************/

.sectionHeader {
    font-size: 1.3em;
    font-weight: normal;
    color: #0206AE;
    background-color: #d9d9d9;
    padding-left: 10px;
    line-height: 34px;
    letter-spacing: 1px;
    margin: 0;
}
.colHeader {
    color: #0206AE;
    background-color: #e7e7e7;
    font-size: 11px;
    line-height: 20px;
    font-weight: normal;
    padding-left: 10px;
    letter-spacing: 0.5px;
    margin: 2px 0 0 0;
}
.navigator {
    color: #bcbcbc;
    background-color: #444444;
    padding-left: 10px;
    font-weight: bold;
    text-align: right;
    line-height: 18px;
}
    
/* /////////////// Article Manager /////////////// */

.forApproval {
    font-weight: bold;
    color: #ff0000;
}
.approved {
    font-weight: bold;
    color: #ff9933;
}
.published {
    font-weight: bold;
    color: #00cc00;
}
.archived {
    font-weight: bold;
    color:  #909090;
}  

/******************************************************************************/
/*                                  LEGACY CSS                                */
/*                                                                            */
/* Note: I am removing elements from here as I replace them with new CSS      */
/*       elements above.  Eventually, there shouldn't be any CSS left here.   */
/******************************************************************************/

/* /////////////// Table modifiers  /////////////// */

.fieldHeader {
    background-image: url(bg_header.png);
    background-repeat: repeat-x;
    color: #000000;
}
.fieldHeader2 {
    background-image: url(bg_header2.png);
    background-repeat: repeat-x;
    color: #000000;
    font-size: 12px;
}
.fieldHeader span, .fieldHeader2 span {
    font-size: 18px;
    padding: 0 6px 0 0;
}
.fieldName, .fieldNameWrap {
    background-color: #e7e7e7;
    color: #0206AE;
    font-weight: bold;
    text-align: left;
    width: 150px;
    border: 1px outset;
}
.fieldName {
    white-space: nowrap;
}
.fieldValue {
    background-color: #F4F7FA;
    color: #000000;
    font: 11px Verdana, Arial, Helvetica, sans-serif;
    line-height: 16px;
}
.fieldValue2 {
    background-color: #F8FDF8;
    color: #0206AE;
    font: 11px Verdana, Arial, Helvetica, sans-serif;
    line-height: 16px;
}
.fieldValue a, .fieldValue a:visited, .fieldValue2 a, .fieldValue2 a:visited {
    font-weight: normal;
}
.subfieldName {
    background-color: #F5F5F5;
    color: #0206AE;
    font-weight: bold;
    text-align: center;
    border: 1px outset;
}
.subfieldValue {
    background-color: #F4F7FA;
    color: #000000;
    font: 11px Verdana, Arial, Helvetica, sans-serif;
    line-height: 16px;
}
.subfieldValue a, .subfieldValue a:visited {
    font-weight: normal;
}
.fieldNote {
    font-weight: bold;
    font-size: 12px;
    white-space: nowrap;
}
.Search {
    background-image: url(bg_search.png);
    background-repeat: repeat-x;
    color: #0206AE;
    text-align: center;
}
.List {
    background-image: url(bg_list.png);
    background-repeat: repeat-x;
    text-align: center;
    font-size: 12px;
}
.newsItem {
    margin: 0 10%;
    border: 2px dotted #0206AE;
    padding: 10px 10px;
    background-color: #ffffcc;
    color: #0206AE;
}
.newsTitle {
    font-weight: bold;
    font-size: 16px;
    text-align: center;
    padding-bottom: 10px;
}
.newsBody {
    font-size: 11px;
}
.results {
    background-color: #F5F5F5;
    background-image: url(bg_pager.png);
    background-repeat: repeat-x;
    color: #9dcdfe;
    font-weight: bold;
    font-size: 12px;
    border: 2px outset;
    padding: 1px 40px;
    white-space: nowrap;
}
.mandatory {
    background-color: #F5F5F5;
    background-image: url(bg_mandatory.png);
    background-repeat: repeat-x;
    font-weight: bold;
    font-size: 10px;
    text-align: center;
    white-space: nowrap;
    border: 1px outset;
}
fieldset {
    color: #0206AE;
    font-size: 1.1em;
    font-weight: bold;
}
legend {
    color: #0206AE;
}

/* /////////////// Links  /////////////// */

.linkCrumbsAlt1 {
    text-decoration: none;
    color: #0206AE;
    font-weight: normal;
    letter-spacing: 0.5px;
}
.linkCrumbsAlt1:hover {
    text-decoration: underline;
    color: #0206AE;
}

/* /////////////// Various /////////////// */

div.pinstripe table {
    background-color: #e7e7e7;
    width: 90%;
}
div.pinstripe td {
    background-color: #ffffff;
}
div.pinstripe img {
    padding: 10px;
}
div.pinstripe button {
    padding: 10px 0;
}
.noBorder {
    border: none;
    font-size: 10px;
}
ul.noindent {
    list-style-image: url('http://www2.ba.itb.cnr.it/UTRSite/themes/default/images/bullet.gif');
    margin: 5px 0 5px 10px;
    padding-left: 5px;
}
ul.bullets li {
    list-style-image: url('http://www2.ba.itb.cnr.it/UTRSite/themes/default/images/bullet.gif');
    margin-left: -15px; 
}
div.bullets {
    margin: -15px 0;
}
.pager {
    background-color: #0206AE;
    background-image: url(bg_pager.png);
    background-repeat: repeat-x;
    color: #9dcdfe;
    white-space: nowrap; 
    width: 70%; 
    margin: 0 auto 5px auto; 
    border: 2px outset;
}
.pager a, .pager a:visited {
    color: #9dcdfe;
    text-decoration: none;
    padding: 1px 4px;
}
.pager a:hover {
    background-color: #0206AE;
    background-image: url(bg_pager2.png);
    color: #FFFFFF;
    text-decoration: none;
    padding: 1px 4px;
}
.pager u {
    background-color: #0206AE;
    background-image: url(bg_pager2.png);
    color: #0206AE;
    text-decoration: none;
    padding: 1px 4px;
}
.narrowButton {
    text-align: center;
    width: 100px;
    line-height: 14px;
    letter-spacing: 1px;
    font-weight: bold;
}
.wideButton {
    text-align: center;
    width: 140px;
    line-height: 14px;
    letter-spacing: 1px;
    font-weight: bold;
}
.linkButton {
    color: #0206AE;
    text-align: center;
    width: 80px;
    line-height: 14px;
    font-size: 11px;
    font-weight: bold;
}
.link2Button {
    color: #0206AE;
    text-align: center;
    width: 50px;
    line-height: 14px;
    font-size: 11px;
    font-weight: bold;
}
.moreButton {
    color: #0206AE;
    text-align: center;
    width: 100px;
    line-height: 14px;
    letter-spacing: 1px;
    font-weight: bold;
}
.addButton {
    color: #008000;
    text-align: center;
    width: 120px;
    font-weight: bold;
}
.deleteButton {
    color: #C00000;
    text-align: center;
    width: 80px;
    font-weight: bold;
}
.searchButton {
    color: #0206AE;
    text-align: center;
    width: 120px;
    font-weight: bold;
}
.message {
    margin: 0 auto;
    border: 2px dotted #ff9600;
    background-color: #ffff99;
    text-align: center;
    width: 50%;
    font-weight: bold;
}
.messageContainer {
    margin: 0 auto 5px auto;
    width: 300px;
}
.messageErrorTitle {
    background-color: #C00000;
    background-image: url(bg_error.png);
    background-repeat: repeat-x;
    color: #ffffff;
    font-weight: bold;
    letter-spacing: 1px;
    text-align: center;
    text-transform: uppercase;
}
.messageError {
    border: 2px dotted #C00000;
    border-top: 1px solid #C000000;
    color: #999999;
    background-color: #ffff99;
    text-align: left;
    padding: 0 0 0 10px;
}
.bgnd {
    background-color: #e5f1ff;
    border: 1px solid #aaaaaa;
}
.treeMenuDefault {
    font-size: 11px;
}

/* /////////////// Tooltips /////////////// */

span.tipOwner {
    position: relative; 
    cursor: hand;
    text-decoration: none;
    float: right;
}
span.tipOwner span.tipText {
    display: none;
    position: absolute;
    top: 0;
    padding: 5px;
    border: 1px solid #04A3EC;
    background-color: #e7e7e7;
    color: #000000;
    font-size: 100%;
    text-align: center;
    width: 13em;
    -moz-opacity: 0.85;
    filter: alpha(opacity=85);
    filter: progid: DXImageTransform.Microsoft.Alpha(opacity=85);
}
/* show/hide: for IE via .htc file, for non-IE via :hover pseudo */
span.tipOwner:hover span.tipText {
    display: block;
}
span.tipOwner {
//  behavior: url(themes/default/css/tooltipHover.htc);
    behavior: url(http://www2.ba.itb.cnr.it/UTRSite/themes/default/css/tooltipHover.htc);
}

