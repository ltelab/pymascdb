.. _examples:

Examples
=======================================
Examples here TODO

.. raw:: html

    <embed>
	<style type="text/css">
.tg  {border-collapse:collapse;border-spacing:0;}
.tg td{border-color:black;border-style:solid;border-width:1px;font-family:Arial, sans-serif;font-size:14px;
  overflow:hidden;padding:10px 5px;word-break:normal;}
.tg th{border-color:black;border-style:solid;border-width:1px;font-family:Arial, sans-serif;font-size:14px;
  font-weight:normal;overflow:hidden;padding:10px 5px;word-break:normal;}
.tg .tg-c3ow{border-color:inherit;text-align:center;vertical-align:top}
.tg .tg-0pky{border-color:inherit;text-align:left;vertical-align:top}
</style>
<table class="tg">
<thead>
  <tr>
    <th class="tg-c3ow">Parameter</th>
    <th class="tg-0pky">Units</th>
    <th class="tg-0pky">Type</th>
    <th class="tg-0pky">Long name<br></th>
    <th class="tg-0pky">Reference / Format / Algorithm<br></th>
  </tr>
</thead>
<tbody>
  <tr>
    <td class="tg-c3ow">Global information </td>
    <td class="tg-0pky"></td>
    <td class="tg-0pky"></td>
    <td class="tg-0pky"></td>
    <td class="tg-0pky"></td>
  </tr>
  <tr>
    <td class="tg-0pky">datetime</td>
    <td class="tg-0pky"></td>
    <td class="tg-0pky">datetime</td>
    <td class="tg-0pky"></td>
    <td class="tg-0pky"></td>
  </tr>
  <tr>
    <td class="tg-0pky">campaign</td>
    <td class="tg-0pky">-</td>
    <td class="tg-0pky">string</td>
    <td class="tg-0pky">Field campaign string</td>
    <td class="tg-0pky"></td>
  </tr>
  <tr>
    <td class="tg-0pky">latitude</td>
    <td class="tg-0pky">deg</td>
    <td class="tg-0pky">float</td>
    <td class="tg-0pky">WGS84 latitude</td>
    <td class="tg-0pky"></td>
  </tr>
  <tr>
    <td class="tg-0pky">longitude</td>
    <td class="tg-0pky">deg</td>
    <td class="tg-0pky">float</td>
    <td class="tg-0pky">WGS84 longitude</td>
    <td class="tg-0pky"></td>
  </tr>
  <tr>
    <td class="tg-0pky">altitude</td>
    <td class="tg-0pky">m a m s l</td>
    <td class="tg-0pky">float</td>
    <td class="tg-0pky"></td>
    <td class="tg-0pky"></td>
  </tr>
  <tr>
    <td class="tg-c3ow">Flake information</td>
    <td class="tg-0pky"></td>
    <td class="tg-0pky"></td>
    <td class="tg-0pky"></td>
    <td class="tg-0pky"></td>
  </tr>
  <tr>
    <td class="tg-0pky">flake_id</td>
    <td class="tg-0pky">-</td>
    <td class="tg-0pky">string</td>
    <td class="tg-0pky">Unique flake ID</td>
    <td class="tg-0pky">e.g. 2015.02.10_11.55.10_flake_4<br>YYYY.MM.DD_HH.mm.ss_flake_flake_number_tmp</td>
  </tr>
  <tr>
    <td class="tg-0pky">flake_number_tmp</td>
    <td class="tg-0pky">-</td>
    <td class="tg-0pky">string</td>
    <td class="tg-0pky">Temporary flake ID<br>(not unique)</td>
    <td class="tg-0pky"></td>
  </tr>
  <tr>
    <td class="tg-0pky">flake_quality_xhi</td>
    <td class="tg-0pky">-</td>
    <td class="tg-0pky">float</td>
    <td class="tg-0pky">Average quality index<br>(on the three cams)</td>
    <td class="tg-0pky">$\chi$ in Praz et al, 2017 citePraz </td>
  </tr>
  <tr>
    <td class="tg-0pky">flake_fallspeed</td>
    <td class="tg-0pky">m/s</td>
    <td class="tg-0pky">float</td>
    <td class="tg-0pky">Recorded fallspeed</td>
    <td class="tg-0pky"></td>
  </tr>
  <tr>
    <td class="tg-0pky">flake_n_roi</td>
    <td class="tg-0pky">-</td>
    <td class="tg-0pky">int</td>
    <td class="tg-0pky">Average # ROIs<br>(on the three cams)</td>
    <td class="tg-0pky"></td>
  </tr>
  <tr>
    <td class="tg-0pky">flake_Dmax</td>
    <td class="tg-0pky">m</td>
    <td class="tg-0pky">float</td>
    <td class="tg-0pky">Maximum Dmax<br>(on the three cams)</td>
    <td class="tg-0pky">Table A1:4 Praz et al, 2017 citePraz </td>
  </tr>
  <tr>
    <td class="tg-c3ow">Riming estimation information</td>
    <td class="tg-c3ow"></td>
    <td class="tg-c3ow"></td>
    <td class="tg-c3ow"></td>
    <td class="tg-c3ow"></td>
  </tr>
  <tr>
    <td class="tg-0pky">riming_deg_level</td>
    <td class="tg-0pky">-</td>
    <td class="tg-0pky">float</td>
    <td class="tg-0pky">Continuous riming<br>degree level</td>
    <td class="tg-0pky">$R_c$ in Praz et al, 2017 citePraz<br></td>
  </tr>
  <tr>
    <td class="tg-0pky">riming_class_id</td>
    <td class="tg-0pky">-</td>
    <td class="tg-0pky">int</td>
    <td class="tg-0pky">Discrete riming<br>degree class ID</td>
    <td class="tg-0pky">Praz et al, 2017 citePraz<br>0: undefined, 1: unrimed, 2: rimed<br>3: densely-rimed, 4: graupel-like, 5: graupel</td>
  </tr>
  <tr>
    <td class="tg-0pky">riming_class_prob</td>
    <td class="tg-0pky">-</td>
    <td class="tg-0pky">float</td>
    <td class="tg-0pky">Riming classification<br>probability</td>
    <td class="tg-0pky">Praz et al, 2017 citePraz</td>
  </tr>
  <tr>
    <td class="tg-0pky">riming_class_name</td>
    <td class="tg-0pky">-</td>
    <td class="tg-0pky">string</td>
    <td class="tg-0pky">Discrete riming <br>degree class name</td>
    <td class="tg-0pky">See riming_class_id</td>
  </tr>
  <tr>
    <td class="tg-c3ow">Melting estimation information</td>
    <td class="tg-0pky"></td>
    <td class="tg-0pky"></td>
    <td class="tg-0pky"></td>
    <td class="tg-0pky"></td>
  </tr>
  <tr>
    <td class="tg-0pky">melting_class_id</td>
    <td class="tg-0pky">-</td>
    <td class="tg-0pky">int</td>
    <td class="tg-0pky">Discrete melting <br>class ID</td>
    <td class="tg-0pky">Praz et al, 2017 citePraz<br>0: dry, 1: melting</td>
  </tr>
  <tr>
    <td class="tg-0pky">melting_prob</td>
    <td class="tg-0pky">-</td>
    <td class="tg-0pky">float</td>
    <td class="tg-0pky">Melting probability </td>
    <td class="tg-0pky">Praz et al, 2017 citePraz<br>If rounded, it yields melting_class_id</td>
  </tr>
  <tr>
    <td class="tg-0pky">melting_class_name</td>
    <td class="tg-0pky">-</td>
    <td class="tg-0pky">string</td>
    <td class="tg-0pky">Discrete melting<br>class name</td>
    <td class="tg-0pky">See melting_class_id</td>
  </tr>
  <tr>
    <td class="tg-c3ow">Hydrometeor type estimation information</td>
    <td class="tg-0pky"></td>
    <td class="tg-0pky"></td>
    <td class="tg-0pky"></td>
    <td class="tg-0pky"></td>
  </tr>
  <tr>
    <td class="tg-0pky">snowflake_class_name</td>
    <td class="tg-0pky">-</td>
    <td class="tg-0pky">string</td>
    <td class="tg-0pky">Hydrometeor <br>class name</td>
    <td class="tg-0pky">Praz et al, 2017 citePraz<br>See snowflake_class_id<br></td>
  </tr>
  <tr>
    <td class="tg-0pky">snowflake_class_id</td>
    <td class="tg-0pky">-</td>
    <td class="tg-0pky">int</td>
    <td class="tg-0pky">Hydrometeor<br>class ID</td>
    <td class="tg-0pky">Praz et al, 2017 citePraz<br>1: small_particle, 2: columnar_crystal, <br>3: planar_crystal, 4: aggregate,<br>5: graupel, 6: columnar_planar_combination</td>
  </tr>
  <tr>
    <td class="tg-0pky">snowflake_class_prob</td>
    <td class="tg-0pky">-</td>
    <td class="tg-0pky">float</td>
    <td class="tg-0pky">Classification <br>probability</td>
    <td class="tg-0pky"></td>
  </tr>
  <tr>
    <td class="tg-c3ow">3D reconstruction / mass estimation</td>
    <td class="tg-0pky"></td>
    <td class="tg-0pky"></td>
    <td class="tg-0pky"></td>
    <td class="tg-0pky"></td>
  </tr>
  <tr>
    <td class="tg-0pky">gan3d_mass</td>
    <td class="tg-0pky">kg</td>
    <td class="tg-0pky">float</td>
    <td class="tg-0pky">Estimated mass</td>
    <td class="tg-0pky">Leinonen et al, 2021 citeLeinonen</td>
  </tr>
  <tr>
    <td class="tg-0pky">gan3d_volume</td>
    <td class="tg-0pky">m$^3$</td>
    <td class="tg-0pky">float</td>
    <td class="tg-0pky">Estimated volume</td>
    <td class="tg-0pky">Leinonen et al, 2021 citeLeinonen</td>
  </tr>
  <tr>
    <td class="tg-0pky">gan3d_gyration</td>
    <td class="tg-0pky">m</td>
    <td class="tg-0pky">float</td>
    <td class="tg-0pky">Estimated gyration<br>radius</td>
    <td class="tg-0pky">Leinonen et al, 2021 citeLeinonen</td>
  </tr>
  <tr>
    <td class="tg-c3ow">Co-located environmental information<br></td>
    <td class="tg-0pky"></td>
    <td class="tg-0pky"></td>
    <td class="tg-0pky"></td>
    <td class="tg-0pky"></td>
  </tr>
  <tr>
    <td class="tg-0pky">env_T</td>
    <td class="tg-0pky">deg C</td>
    <td class="tg-0pky">float</td>
    <td class="tg-0pky">Air temperature </td>
    <td class="tg-0pky"></td>
  </tr>
  <tr>
    <td class="tg-0pky">env_P</td>
    <td class="tg-0pky">hPa</td>
    <td class="tg-0pky">float</td>
    <td class="tg-0pky">Pressure</td>
    <td class="tg-0pky"></td>
  </tr>
  <tr>
    <td class="tg-0pky">env_DD</td>
    <td class="tg-0pky">deg</td>
    <td class="tg-0pky">float</td>
    <td class="tg-0pky">Wind direction<br>(North to East)</td>
    <td class="tg-0pky"></td>
  </tr>
  <tr>
    <td class="tg-0pky">env_FF</td>
    <td class="tg-0pky">m/s</td>
    <td class="tg-0pky">float</td>
    <td class="tg-0pky">Wind speed</td>
    <td class="tg-0pky"></td>
  </tr>
  <tr>
    <td class="tg-0pky">env_RH</td>
    <td class="tg-0pky">\%</td>
    <td class="tg-0pky">float</td>
    <td class="tg-0pky">Relative humidity</td>
    <td class="tg-0pky"></td>
  </tr>
  <tr>
    <td class="tg-c3ow">Blowing snow estimation</td>
    <td class="tg-0pky"></td>
    <td class="tg-0pky"></td>
    <td class="tg-0pky"></td>
    <td class="tg-0pky"></td>
  </tr>
  <tr>
    <td class="tg-0pky">bs_normalized_angle</td>
    <td class="tg-0pky">-</td>
    <td class="tg-0pky">float</td>
    <td class="tg-0pky">Blowing Snow<br>normalized angle</td>
    <td class="tg-0pky">Schaer et al 2020, citeSchaer<br>Pure precip. if $&lt; 0.193$, Pure BS if $&gt; 0.881$</td>
  </tr>
  <tr>
    <td class="tg-0pky">bs_mixing_ind</td>
    <td class="tg-0pky">-</td>
    <td class="tg-0pky">float</td>
    <td class="tg-0pky">Blowing snow<br>mixing index</td>
    <td class="tg-0pky">Schaer et al 2020, citeSchaer<br>Only defined in mixed BS/precip environments</td>
  </tr>
  <tr>
    <td class="tg-0pky">bs_precip_class_name</td>
    <td class="tg-0pky">-</td>
    <td class="tg-0pky">string</td>
    <td class="tg-0pky">Blowing snow <br>class name</td>
    <td class="tg-0pky">Schaer et al 2020, citeSchaer<br>See bs_precip_class_id</td>
  </tr>
  <tr>
    <td class="tg-0pky">bs_precip_class_id</td>
    <td class="tg-0pky">-</td>
    <td class="tg-0pky">int</td>
    <td class="tg-0pky">Blowing snow<br>class ID</td>
    <td class="tg-0pky">Schaer et al 2020, citeSchaer<br>0: undefined, 1: precip, 2: mixed,<br>3: blowing_snow<br></td>
  </tr>
</tbody>
</table>
      
    </embed>
