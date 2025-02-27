#!/bin/sh

rm -rf ./www2.ba.itb.cnr.it
rm -rf ./UTRSite
wget -r http://www2.ba.itb.cnr.it/UTRSite/index.php/
cp -rf ./www2.ba.itb.cnr.it/UTRSite/index.php/UTRSite\ signal/Signal/action/view/frmUTRSiteID ./UTRSite
