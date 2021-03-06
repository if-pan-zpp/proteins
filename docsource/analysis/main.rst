Main/global code
================

Overview
--------
The code (``code_notes.f``, which is ``cg.f`` but linted and with added inline comments), takes a form of an initial "imperative" header, i.e. global/main code, followed by a sequence of subroutine definitions. For the sake of focus we shall extract the main part into ``splices/main.f``. It takes ~2.5k lines, and some general comments are provided.

Preparation
-----------

First, variables are declared (lines 1-72), and the initial values are specified (lines 74-273). Then, an input file is read: the program itself is executed, according to README.txt, as ``./cg inputfile`` (line 37), so the ``inputfile`` string is extracted (lines 275-279). Apparently the switch at line 275 indicates that we may choose not specify the input file. At any rate, lines 278-621 parse the key-value pairs in a standard fashion, inexcept for:

- setting ``mapfile``, which also sets ``lwritemap`` to true;
- setting ``cmapf``, which also sets ``lcmap`` to true;
- setting ``lecperm``, which also initializes? ``coul`` to 210;
- setting ``displace`` -- this one, as described in README.txt, takes the numbers of the chains (``iprota``, ``iprotb``), and the desired displacement distance ``away``; also, ``ldisp`` is set to true;
- setting ``paramfile``, which also sets ``lparam`` to true;
- setting ``period`` or ``omega``, which are physically the same (except for scaling etc.);
- setting ``bbrm``, ``ssrm``, ``bsrm``, ``i4rm``. (TODO: Write a separate section on ``sigma1``);
- setting ``temp``, which sets to its value variables ``tstart`` and ``tend`` -- there is no variable ``temp`` on its own;
- setting ``file``, which sets ``outfile``, ``mapfile`` and ``savfile`` to the name but with suffixes.

Then (line 622) variables are loaded from a parameter file (in a separate subroutine ``load_paramfile``, why aren't the input file variables also in a separate subroutine?).

Angles ``the0`` and ``phi0`` are reset, and ``lfrompdb`` is cleared (lines 624-626). This variable, I presume, is used for computing angle potentials, and ``lsimpang=.true.`` should mark the use of simple angle potentials. If ``lpbc`` is true, periodic boundary conditions are used on all planes (is this "turn on/off for each axis" functionality ever used?) (lines 628-632). Lines 634-641 are for the most part documented. Then, for some reason ``#I,I+2 CONTACTS PURELY REPULSIVE`` is written (unless there's some jump there? I don't see it). The relevant variables are scaled (in the code, according to CPC14.pdf, units are 5 A for distance, :math`\varepsilon \approxeq` 1.5 kcal/mol for energy, \approxeq 1ns for time are used).

Then, we "set up the table of temperature" in lines 689-696. More specifically, we can specify (via the use of different ``tstart`` and ``tend`` in the input file) to have the ambient temperature of the solution change, linearly in at most 150 stages (``tstep`` declares the desired temperature change in one step). In this sense, the temperature progression is not linear per se, but rather "stratified", like a ceil function.

The next section is titled as such (line 701). The lines 702-719 about loading from either PDB file or a sequence file. In the former case, ``load_protein`` is invoked; in the latter, ``load_sequence`` followed by ``confstart``. Also, ``x0..`` and ``xn..`` are made to be the same. Let us speak of these, because it's confusing. ``x0..`` and ``xn..`` are supposed to represent, respectively, *actual* and *native* coordinates (line 942). The question, then, is why:

- these are, in some sense, reset at one point (lines 1254-1267);
- ``xn..`` are used for computing force when using AFM protocol (lines 1347-1349).

In general, the matter of the stages before the simulation proper is strange (in my estimation); we will return to the matter later.

Having the sequences and coordinates, we construct contact maps (lines 721-736). Lines 737-740 set ``ipwn``, which is supposed to be the number of residues "sticked to wall" (whatever that means) (line 229). As for the formula (line 738), I don't really understand why it's the way it is (especially with the cubic roots and whatnot). Then (lines 741) we build "titin", but from what I have seen this part is essentially not used. In lines 749-758, residue masses are set. Next part (lines 759-792) set up some auxiliary variables, mostly regarding the "simulation box". Line 793 sets appropriate contact type for the setting where "residues are connected to walls one by one"; I will for now skip anything regarding the walls. Then (line 794) Verlet list(s) are updated.

In lines 795-809 there's something about "contact orders". ``krist`` is an array containing repulsive contacts (presumably derived from the distances, the only place where ``krist`` is set are lines 4462-4465 of ``code_notes.f``, in ``update_verlet_list``). From the analysis of these lines, ``corder`` as of line 806 is the relative distance in residue indices between repulsive contacts. Interestingly, if ``icor`` is 0, a message ``#NO NON-REPULSIVE CONTACTS`` is written, but in ``update_verlet_list`` where new entries are added, there's ``! REPULSIVE CONTACT``, so I don't know. For that matter, I haven't seen ``corder`` used in a "functionally relevant way", only printed, for example in line 2137. Lines 810-821 are much the same, but about "native contacts".

Lines 823-854 are about disulfide bridges: for every one, residues are checked whether they are indeed cysteines, and if these bonds are not to be treated as normal (``lsslj``) or exclusive (``llsselj``) L-J contacts, the relevant contact (from ``klist``) is set to have a "type" +=631, for some reason.

Lines 857-878 do some logging, it's sort of irrelevant. Lines 882-894 do some further scaling (why it's not in the previous scaling section is unclear, perhaps the subroutines set these values, but then why aren't the values set by subroutines made scaled by default?).

Chirality potentials are optionally computed (line 901), much like the native angles (if the residues were loaded from a PDB file), in line 916. In lines 918-938, according to the message, we "disable native contacts (I, I+2)". In CPC14.pdf (p. 5?), interactions between so close residues are regarded as "local", as opposed to "non-local" interactions which the contact maps service - thus we remove such contacts, if they occur. This (lines 924-929) is also the place where dynamic SS bonds are logged (note that ``ijdiff`` is set to zero, so unless ``.not.lconect(i)``, the later code is ignored).

Lines 940-957 do some further initialization. For some reason, glutamines are counted; I checked where ``ngln`` may be used, and it's apparently just logged for whatever (presumably debug) reason. Lines 958-971 set up the sidechain types if the chains are to be charged; further elaboration may be required, but I don't know where precisely the electrostatic forces are used so whatever. Lines 973-985 chiefly set up ``b`` array (which makes it hard to search for references); this array contains distances of consecutive Ca-Ca bonds.

Then (lines 987-1012), we essentially prepare other stuff (again, why was the "free-floating" stuff also not extracted into subroutines? Not that merely extracting code and putting it in a subroutine is good, but for example ``gopotential`` is like that, so it's not like it's consistent). In order, we call ``gopotential`` (and potentially print it), ``prepare``, ``evalgo``, ``evalwall``, ``evalimproper`` (or ``evalcpot`` depending on ``lpid``), ``evalangles``, ``eval_chirality`` (not sure what line 1008 is doing). Lines 1013-1070 are preoccupied with logging a lot of stuff. Then we rescale another parameters (another scaling section) (lines 1073-1097). Lines 1102-1111 print out the successive temperatures (see ``tstart``, ``tend``).

Main loop
---------
Lines 1113-1151 set some parameters up. Then (starts at line 1156), we get to the "meat" of the program, i.e. looping over the ``ntraj`` starting configurations. Lines 1157-1236 for the most part "merely" reset the variables (presumably across the various trajectories); some assignments are annotated. (Note: an explanation for the role of ``*resist`` variables is sorely needed).

Then, we have a "starting conformation" phase (lines 1238-1267). Here we have (I think first?) example of simulation protocols. If we perform a folding simulation (``lconftm`` or ``lmedian`` true), then we want the "straight line" starting configuration. For ``lwarmup`` true, we set the positions to the native ones *for the first trajectory*, and otherwise for every trajectory if the residues have a native structure. If not, we call ``confstart`` (which initializes the locations), and save them into the native coordinates array ``xn..``.

In the next order, we (optionally) displace the chains (line 1269). Lines 1271-1300 are about the variables for periodic boundary conditions (hard to say at this moment how they are used). Lines 1301-1338 compute wall stuff. We apparently set up arrays for sorting residues with respect to z-axis (lines 1305-1311); if I remember correctly, it's so that we can sort the residues for attaching checking whether they should be attached to xy-plane walls *somewhere*. Lines 1312-1337 are triggered if we want to simulate squeezing the simulation box to reach desired density. Lines 1339-1344 save original simulation box planes, which are used later. Then (1346-1357) "the direction of force" is computed (for AFM stretching simulation). Speaking of AFM stretching simulation, we can set ``lvelo`` from the input file, but not ``lforce``, seems like an unfinished feature. If said stretching is to have constant velocity (``lvelo``), we set up max residue length (of course, if the residue is already maximally extended, there's no point in stretching), and ``xpul..`` variables (places where the "tip" is attached?).

Then (line 1375-) we initialize the velocities, set ``kb0`` to zero (not sure what this variable is for exactly, as with ``sep0``). In lines 1376-1397 we "zero the connection tables"; from the commented-out line 956, and the later occurences of the variable, ``khbful`` seems to have something to do with the current number of bonds that a residue may have with other residues. ``nei`` is the neighbour counter (see line 952). If dynamic SS bonds are disabled, we set them up from the array of static ones in lines 1388-1397.

Lines 1398-1460 restore the state from a reset file; other than the fact that it's done manually, there's (in my estimation) nothing sophisticated about it. Lines 1462-1486 in some sense replicate lines 997-1009, so the question would be why were these invoked before the trajectory loop (perhaps to initialize some auxiliary variables? not sure). Then we scale accelerations for some reason (lines 1489-1493); in lines 1496-1535 just log stuff. In lines 1537-1568, we simulate pulling with an AFM for some time (``mpull``) before finally stopping. In some sense this is interesting, because such pulling shouldn't necessarily be any different from normal AFM simulation (except for being limited), and because the code is fairly concise -- we may therefrom gain more insight into the role each part plays in the "main simulation" part below.

Thus we enter "main loop of simulation". Lines 1574-1576 indicate, that ``kb0`` is the initial "frame number" for the simulation proper, and that ``kb`` is the current frame. Lines 1578-1612 follow a previously-seen formula (like 1541-1559), but now we account for ``lmass`` option, and there are modifications in 1601-1604. More specifically, ``lsldh`` is for "slowly turning on dihedral potential.

Lines 1613-1903 are for "equilibration". Lines 1613-1872 handle the equilibration frames; for the most part this phase seems to apply only if the simulation protocol includes walls (lines 1615-1849), followed by lines 1851-1854 which I'm not sure what they do, followed by "applying" force in lines 1856-1871 (for AFM protocol, from what I see). (TODO: finish this part, as I can't be bothered to analyze it now and it isn't relevant for the first prototype except for the lines 1851-1854 perhaps).

Lines 1905-1933 are supposed to "measure the force to walls (A and D)", which in particular has a lot of stuff about ``*resist`` variables, so frankly I don't know what is going on here. It may, however, be helpful in understanding these variables. Lines 1935-1938 handle counting the formation of native contacts; in particular, if all are formed, we break the simulation (this might be slightly modified depending on the type of SS formation). Line 1939 skips certain simulation steps (interestingly it's only *after* the whole equilibration thing). Lines 1941-1943 seem unfinished. Lines 1945-1971 are annotated by default, wherefore I shall omit their description here. Lines 1973-1992 handle that custom trick that Mioduszewski spoke of for reducing the number of Verlet list recomputations needed (?).

Lines 1996-2010 refer to "recording average force when pulling with AFM with constant velocity" (and line 2001, initializing ``ree``, may be of use for stopping simulation if the chain is maximally stretched, as mentioned before). Lines 2012-2027 apparently do something similar? but for the walls, I guess, I don't know. Lines 2029-2201 are apparently "only" for logging stuff to output files, so I'd say they aren't really relevant. Lines 2203-2209 "compute thermodynamic properties" (if the docs for ``lthermo`` are anything to go by); lines 2211-2213 have annotations, and 2214-2218 are I believe about recording median folding times (in particular ``goto 533`` returns to the beginning of the simulation loop proper at line 1575). The rest is logging stuff and loop structures.

Remarks
-------
- Given the very end of the "main" section, specifically the line 2420 with ``2000 continue`` (which refers to the temperature loop), and the comment above, it would "seem" that in general the loop structure is over *temperatures*, then over trajectories, and then over simulation steps. The fact that it's over *temperatures* is somewhat strange in my estimation, so perhaps I am incorrect in this analysis.
