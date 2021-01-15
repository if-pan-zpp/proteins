

## Notes:



*   lpdb = true restricts usage to structured proteins.


## Special interactions


<table>
  <tr>
   <td>Lp
   </td>
   <td>Physical process
   </td>
   <td>Ref
   </td>
   <td>Ref in code
   </td>
   <td>Conditions
   </td>
   <td>Comments
   </td>
  </tr>
  <tr>
   <td>1
   </td>
   <td>Damping
   </td>
   <td>2.
   </td>
   <td>lang, lang_mass
   </td>
   <td>-
   </td>
   <td>
   </td>
  </tr>
  <tr>
   <td>2
   </td>
   <td>Thermal noise
   </td>
   <td>2.
   </td>
   <td>lang, lang_mass
   </td>
   <td>-
   </td>
   <td>
   </td>
  </tr>
</table>



## Local interactions


<table>
  <tr>
   <td>Lp
   </td>
   <td>Physical process
   </td>
   <td>Ref
   </td>
   <td>Ref in code
   </td>
   <td>Conditions
   </td>
   <td>Comments
   </td>
  </tr>
  <tr>
   <td>3
   </td>
   <td>Harmonic tethering between (i, i+1) within one chain.
   </td>
   <td>2. 
   </td>
   <td>evalgo
   </td>
   <td>-
   </td>
   <td>
   </td>
  </tr>
  <tr>
   <td>4
   </td>
   <td>Chirality potential
   </td>
   <td>3.1
   </td>
   <td>eval_chirality
   </td>
   <td>lpdb=T
<p>
lchiral=T
   </td>
   <td>
   </td>
  </tr>
  <tr>
   <td>5
   </td>
   <td>Bond angle (structured)
   </td>
   <td>3.2.1
   </td>
   <td>evalangles
   </td>
   <td>langle=T
<p>
Locally all lfrompdb(i)=T
   </td>
   <td>TODO: describe lfrompdb exactly
   </td>
  </tr>
  <tr>
   <td>6
   </td>
   <td>Bond angle (unstructured, tabularized)
   </td>
   <td>3.2.3
   </td>
   <td>evalangles
   </td>
   <td>langle=T \
Locally some lfrompdb(i)=F
<p>
lenetab=T
   </td>
   <td>
   </td>
  </tr>
  <tr>
   <td>7
   </td>
   <td>Bond angle (unstructured, statistical)
   </td>
   <td>3.2.2
   </td>
   <td>evalangles
   </td>
   <td>langle=T
<p>
Locally some lfrompdb(i)=F
<p>
lenetab=F
   </td>
   <td>It’s called statistical, because its coefficients are obtained from fitting the formulas to potentials resulting from inverse Boltzmann method applied to random coil database.
   </td>
  </tr>
  <tr>
   <td>8
   </td>
   <td>Dihedral angle (structured, harmonic)
   </td>
   <td>3.2.1
   </td>
   <td>evalangles
   </td>
   <td>langle=T
<p>
ldi=T
<p>
Locally all lfrompdb(i)=T
<p>
ldisimp=T
   </td>
   <td>‘ldi’ is the same as ‘ldih’, but ‘ldi’ is used in code and ‘ldih’ in documentation.
   </td>
  </tr>
  <tr>
   <td>9
   </td>
   <td>Dihedral angle (structured, non-harmonic)
   </td>
   <td>3.2.1
   </td>
   <td>evalangles
   </td>
   <td>langle=T
<p>
ldi=T
<p>
Locally all lfrompdb(i)=T
<p>
ldisimp=F
   </td>
   <td>
   </td>
  </tr>
  <tr>
   <td>10
   </td>
   <td>Dihedral angle (unstructured, statistical)
   </td>
   <td>3.2.2
   </td>
   <td>evalangles
   </td>
   <td>langle=T
<p>
ldi=T
<p>
Locally some lfrompdb(i)=F
   </td>
   <td>
   </td>
  </tr>
  <tr>
   <td>11
   </td>
   <td>Dihedral angle
<p>
(unstructured, tabularized)
   </td>
   <td>3.2.3
   </td>
   <td>??? missing feature
   </td>
   <td>-
   </td>
   <td>It’s mentioned in README.txt and somehow implied in CPC14/3.2.3, but almost surely isn’t implemented.
   </td>
  </tr>
</table>



## Non-local interactions


<table>
  <tr>
   <td>Lp
   </td>
   <td>Physical process
   </td>
   <td>Ref
   </td>
   <td>Ref in code
   </td>
   <td>Conditions
   </td>
   <td>Comments
   </td>
  </tr>
  <tr>
   <td>12
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
   <td>
   </td>
  </tr>
</table>

