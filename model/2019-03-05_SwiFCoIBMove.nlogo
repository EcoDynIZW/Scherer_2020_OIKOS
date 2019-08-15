;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                                                                                                                                                                                                                           ;;
;;  SwiFCoIBM-Move - Individual-Based Swine Fever Community Model with explicit movement (based on model versions published in Kramer-Schadt et al. 2009 and Lange et al. 2012)                                                                              ;;
;;                                                                                                                                                                                                                                                           ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                                                                                                                                                                                                                           ;;
;;  The model is composed of two major components:                                                                                                                                                                                                           ;;
;;    (1) a wild boar demography model considering seasonal reproduction, herd splitting (dispersal) and mortality and                                                                                                                                       ;;
;;    (2) a Classical Swine Fever (CSF) virus model operating on the emerging wild boar population.                                                                                                                                                          ;;
;;                                                                                                                                                                                                                                                           ;;
;;  Wild boar population density and structure are affected by the disease via virus-induced mortality and litter size depression.                                                                                                                           ;;
;;  The crucial model entity is the wil boar individual, characterised by age in weeks (one week represents the approx. CSF incubation time; Artois et al. 2002, Moenning et al. 2003),                                                                      ;;
;;    resulting in age classes piglets (age <= 8 months), yearlings (subadults; 8 months < age <= 2 years) and adults (age > 2 years).                                                                                                                       ;;
;;  Each host has a location, which denotes its home range cell and the individual's family group (
;;
;;                                                                                                                                                                                                                                                           ;;
;;  Movement: Adult female movement is restricted to the 'activity zone' which equals their home range cell (staying strategy). Movement is similar for individuals of the sounder (social unit organised around adult females and their offspring).         ;;
;;            Subadult females split in specified weeks of the year from their herds searching for an unoccupied home range cell (see Procedure "HerdSplit").                                                                                                ;;
;;            Subadult males disperse from natal family groups in subgroups at an age the males begin to seek mates (~5 months), joining other groups in the surrounding area.                                                                               ;;
;;            Adult and elderly males tend to range solitary outside of their activity zone (ranging strategy).                                                                                                                                              ;;
;;                                                                                                                                                                                                                                                           ;;
;;  Landscape: The model landscape is represented by a grid of 4 km^2 square cells which each encompass a wild boar family group's home range (Leaper et al. 1999).                                                                                          ;;
;;             Each habitat cell is characterised by a breeding capacity (quality), denoting habitat quality by the number of female boars that are allowed to breed, representing density regulation in the model.                                          ;;
;;                                                                                                                                                                                                                                                           ;;
;;                                                                                                                                                                                                                                                           ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                                                                                                                                                                                                                           ;;
;;  STATE VARIABLES                                                                                                                                                                                                                                          ;;
;;                                                                                                                                                                                                                                                           ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


; NAME                            TYPE/RANGE     DEFAULT        DESCRIPTION

globals
[
  habitat                     ;   [set]          -              habitat cells (cells with quality > 0)
  matrix                      ;   [set]          -              matrix cells (cells with quality = 0)
  ;mean_quality               ;   [#]            9              maximum value for habitat quality to define range of habitat quality (quality between 0 and mean_quality)
  ;year_release               ;   [#]            6              year in which an randomly individual becomes infected (range of week_release)
  week_release                ;   [#]            261-312        random week in year_release of virus release following an expontential distribution with lambda = 1.5
  ;herd_prop                  ;   [%]            100            proportion of habitat patches occupied by herds
  ;release_fct                ;   [#, >=1]       3              multiplies breeding capacity of each occupied patch to determine number of boars
  max_age                     ;   [#]            572 weeks      maximum age in weeks
  ;longevity                  ;   [#]            11 years       maximum age in years (max_age / 52)
  ;age_blur                   ;   [#]            6 weeks        stochastic deviation from age class thresholds in weeks (in both directions (age class threshold +/- age_blur) following a normal distribution
  ;dist_disp                  ;   [#]            3 cells        dispersal distance for subadults
  survival_prob               ;   [%]            0.65           annual survival rate of adults and subadults (mean = 0.65, min = 0.4)
  survival_prob_pig           ;   [%]            0.5            annual survival rate of piglets (mean = 0.5, min = 0.1)
  ;case_fatality              ;   [%]            0.8            general risk for getting lethally infected by a disease transmission
  ;b_within                   ;   [#, 0-1]       0.0208         transmission coefficient within the herd
  ;b_between                  ;   [#, 0-1]       0.00208        transmission coefficient between herds
  ;b_move                     ;   [#, 0-1]       -              transmission coefficient during ranging movements
  ;fetal_inf                  ;   [#, 0-1]       0.5            probability of pre-natal infection
  ;fert_red                   ;   [#, 0-1]       0.625          reduction in fertility of infected females
  ;mue                        ;   [#]            3 weeks        mean of exponential survival time distribution for lethally infected individuals
  ;mue_max                    ;   [#]            10             factor for truncating maximum infectious period (mue * mue_max)
  tick_leth_max               ;   [#]            based on mue   maximum step lethally infected individuals survive (based on mue * mue_max)
  ;t_trans                    ;   [#]            1 week         infectious period of transiently infected individuals
  ;t_anti                     ;   [#]            15 weeks       maximum persistence of maternal antibodies (Depner et al. 2000)
  repro_prob                  ;   [#]            list           monthly reproduction probability
  week                        ;   [#, 1-52]      0              current week of the year
  ;run_years                  ;   [#]            100            number of simulated years during one run
  ;fem_prob                   ;   [%]            0.5            sex ratio (0: only males; 1: only females)
  new_leth                    ;   [#]            0              number of individuals which got lethally infected during this timestep
  new_trans                   ;   [#]            0              number of individuals which got transient infected during this timestep
  ;q                          ;   [#, 0-1]       -              probability to move straight to target cells
  list_n_max                  ;   [list]         []             input data: sums of the number of patches adjacent to each patch in a circle of the same area as the infected/infectious area (see Report_Fragmentation)
  list_repro                  ;   [list]         []             input data: cumulative reproduction probabilities for each month
  rep                         ;   [#]            -              number of repetitions in case of running BehaviourSpace experiments
  ;seed                       ;   [#]            -              seed used if seed_setup = "ON"
  DONE                        ;   [boolean]      0              1 if stop conditions are TRUE

  ;; OUTPUT
  runtime                     ;   [min]          0              counting minutes since starting setup
  dist_inf                    ;   [cells]        0              maximum spatial distance of disease transmission
  inf_x                       ;   [#]            -              coordinate x of the pathogen-introduced cell
  inf_y                       ;   [#]            -              coordinate y of the pathogen-introduced cell
  week_inf_max                ;   [#]            -              week when infection reached other side of the tubus
  week_last                   ;   [#]            -              week in which the last infected individual was noted
  quality_mean_6              ;   [#]            -              mean habitat quality in a radius of 6 cells from the release cell
  count_init                  ;   [#]            -              number of individuals in release cell
  count_init_roaming          ;   [#]            -              number of adult males in release cell
  F_infected                  ;   [0-1]          -              fragmentation index for infectious cells (patches which contain shedding individuals)
  F_infectious                ;   [0-1]          -              fragmentation index for infected cells (patches which were infected at least one week during the simulation)
  list_pop_struc              ;   [list]         []             current age class distribution of all individuals in years
  list_inf_struc              ;   [list]         []             current age class distribution of infected individuals (epi_stat = esLeth or esTrans) in years
  list_trans_struc            ;   [list]         []             current age class distribution of transient infected individuals (epi_stat = esTrans) in years
  list_leth_struc             ;   [list]         []             current age class distribution of lethally infected individuals (epi_stat = esLeth) in years
  list_mue_curr               ;   [list]         []             current infectious periods of lethally infected individuals (epi_stat = esLeth) in weeks
  list_mue_all                ;   [list]         []             all infectious periods of lethally infected individuals (epi_stat = esLeth) in weeks
  dens_var                    ;   [#]            -              variance of turtles per patch
  dens_roam_var               ;   [#]            -              variance of roaming individuals per patch
  dens_inf_group_var          ;   [#]             -             variance of infected group-living individuals per patch
  dens_inf_roam_var           ;   [#]             -             variance of infected group-roaming individuals per patch
  contacts_med                ;   [#]            0              median number of contacts between infected and susceptible individuals (scaled by visit time = 1 / dist)
  contacts_avg                ;   [#]            0              mean number of contacts between infected and susceptible individuals
  contacts_lwr                ;   [#]            0              25% quartile of contacts between infected and susceptible individuals
  contacts_upr                ;   [#]            0              75% quartile of contacts between infected and susceptible individuals
  contacts_max                ;   [#]            0              maximum number of contacts between infected and susceptible individuals
  contacts_var                ;   [#]            0              variance of contacts between infected and susceptible individuals
  contacts_1                  ;   [#]            0              number of contacts within cells of quality_int = 1
  contacts_2                  ;   [#]            0              number of contacts within cells of quality_int = 2
  contacts_3                  ;   [#]            0              number of contacts within cells of quality_int = 3
  contacts_4                  ;   [#]            0              number of contacts within cells of quality_int = 4
  contacts_5                  ;   [#]            0              number of contacts within cells of quality_int = 5
  contacts_6                  ;   [#]            0              number of contacts within cells of quality_int = 6
  contacts_7                  ;   [#]            0              number of contacts within cells of quality_int = 7
  contacts_8                  ;   [#]            0              number of contacts within cells of quality_int = 8
  contacts_9                  ;   [#]            0              number of contacts within cells of quality_int = 9
  trans_g                     ;   [#]            0              number of within-patch transmissions (caused by b_within)
  trans_w                     ;   [#]            0              number of transmissions (caused by b_move) within the cluster the adult male became infected
  trans_b                     ;   [#]            0              number of transmissions (caused by b_move) in another cluster than the one where the adult male became infected
  patches_0                   ;   [#]            0              number of patches with quality_int = 0 (matrix)
  patches_1                   ;   [#]            0              number of patches with quality_int = 1
  patches_2                   ;   [#]            0              number of patches with quality_int = 2
  patches_3                   ;   [#]            0              number of patches with quality_int = 3
  patches_4                   ;   [#]            0              number of patches with quality_int = 4
  patches_5                   ;   [#]            0              number of patches with quality_int = 5
  patches_6                   ;   [#]            0              number of patches with quality_int = 6
  patches_7                   ;   [#]            0              number of patches with quality_int = 7
  patches_8                   ;   [#]            0              number of patches with quality_int = 8
  patches_9                   ;   [#]            0              number of patches with quality_int = 9
]

turtles-own
[
  age                         ;   [#, 0-max_age] -              individual age in weeks
  is_female                   ;   [boolean]      -              individual sex: female or male
  herd_id                     ;   [#]            -              herd membership
  inf_id                      ;   [#]            -              patch id were the infeciton took place
  family                      ;   [#]            -              [who] of mother
  ind_case_fatality           ;   [%]            ind.           individual mortality rate (case_fatality ^ 2 for adults, sqrt case_fatality for Juveniles)
  age_group                   ;   [string]       -              age group turtle belongs to: piglet, subadult, adult; differs between sexes
  dem_stat                    ;   [string]       dsResident     demographic status: dsResident, dsOffspring, dsDispF, dsRoaming, dsDispM_subadult
  epi_stat                    ;   [string]       esSusc         epidemiological status: esSusc, esTrans, esLeth, esImm, esImmMat
  is_shedding                 ;   [boolean]      0              1 if the individual is infected (epi_stat = esTrans or esLeth
  tick_inf                    ;   [#]            -              week of infection
  tick_leth                   ;   [#]            -              week of death after infection (for epi_stat = "esLeth")
  ind_b_move                  ;   [#, 0-1]       b_move/dist    individual transmission probabilities for ranging individuals based on movement steps of current week (~time spent per cell)
  breed_prob                  ;   [#, 0-1]       -              random number for stochastic effects in monthly reproduction for the given year
  tick_birth                  ;   [#]            -              week of birth
  is_matAB                    ;   [boolean]      0              1 if mother had anti bodies (= was immune)
  is_mother                   ;   [boolean]      0              1 if female had offspring in current year
  is_breeder                  ;   [boolean]      0              1 if female is allowed to breed in current year
  contacts                    ;   []             -              number of contacts of an infected adult males with susceptible individuals
  dist_move                   ;   [#]            -              individual weekly movement distance
  visited                     ;   [#]            -              number of unique visited cells
  tag                         ;   [agent]        -              randomly chosen adult males plotting their movement paths
  tcolor                      ;   [#]            -              time-based color to plot movement paths
]

patches-own
[
  is_occupied                 ;   [boolean]      0              1 if a patch contains a herd of wild boars
  quality                     ;   [#]            [0-9]          habitat quality in terms of breeding capacity = number of females alllowed to breed
  quality_int                 ;   [#]            [0-9]          round-up value of habitat quality (ceiling quality)
  is_habitat                  ;   [boolean]      0              true if quality > 0
  qual_around                 ;   [#]            [0-9]          summed habitat quality of surrounding cells (to decide if movement is possible)
  capacity                    ;   [#]            [0-40]         overall capacity = quality * (mean number of 20 boars per cell / mean quality) = quality * (20 / 4.5) = 4.4445 boars per quality unit
  dens_rel                    ;   [#]            -              relative density of the cell (count turtles-here / capacity)
  count_inf                   ;   [#]            0              number of infectious indiviuals in this habitat cell
  inf_press                   ;   [#]            0              infection pressure based on infected indivudals in and around this patch
  is_infectious               ;   [boolean]      0              1 if patch is infected at the current step
  is_infected                 ;   [boolean]      0              1 if patch has been infected during simulation
  ind_sum                     ;   [#]            0              count of individuals per cell during the first year
  inf_sum                     ;   [#]            0              keeps track of how many times a patch contained one or more infected individuals
  count_repro                 ;   [#]            0              number of reproductive females (subadult + adult) per cell
  count_male                  ;   [#]            0              number of  adult males per cell
  count_disp_fem              ;   [#]            0              number of female dispersers per cell
  count_group                 ;   [#]            0              number of individuals belonging to the group (individuals per cell without adult males)
  disp_cells                  ;   [set]          -              cells in dispersal distance of the focal cell
  W_M                         ;   [#]            -              movement weight (for movement decisions)
  last_visited                ;   [#]            0              to ignore the cell the individual has visited in the previous movement step
  scent                       ;   [#]            -              accumulated visit time of infected adult males
  id                          ;   [#]            0              cluster id of cells
]

extensions [profiler]



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                                                                                                                                                                                                                           ;;
;;  SETUP + GO PROCEDURES                                                                                                                                                                                                                                    ;;
;;                                                                                                                                                                                                                                                           ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;; SETUP ______________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________

to Setup

  ca
  reset-ticks

  set rep 200

  Init_Parameters
  Init_Landscape
  Init_Population

end


;; GO _________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________

to Go

  ;profiler:start

  ;if ticks = 0                [ reset-timer ]

  ifelse ticks mod 52 = 0     [ set week 1 ] [ set week week + 1 ]

  if ticks = week_release - 1 [ Disease_Release ]  ;; - 1 to be synchronous with week of the year

  if ticks >= week_release    [ Disease_Transmission ]

                                Movement_Roaming

  if week = 17                [ Movement_DispMales ]

  if week = 29                [ Movement_DispFemales ]

                                Host_Reproduction

                                Host_Death

                                Host_Ageing

  if ticks >= week_release   [ Disease_Transition ]

                                Update

                                Report_Fragmentation

                                ;set runtime precision (timer / 60) 2

                                tick

  ;profiler:stop
  ;print profiler:report
  ;profiler:reset

end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                                                                                                                                                                                                                           ;;
;;  INITALIZATION PROCEDURES                                                                                                                                                                                                                                 ;;
;;                                                                                                                                                                                                                                                           ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;; SET PARAMETERS _____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
  ;; observer procedure
  ;; Sets all parameters needed for initalisation: Current number of habitat patches and maximum number of habitat patches, thresholds for landscape pattern, case fatality, risk of getting lethally infected, week of virus release.
  ;; Prepare output lists.

to Init_Parameters

  if seed_setup = "ON" [ random-seed seed ]
  if seed_setup = "BS"
  [
    set seed behaviorspace-run-number mod rep + 1
    random-seed seed
  ]

  set max_age longevity * 52
  set week 0
  set week_release ((year_release - 1) * 52) + random 52 + 1   ; year_relase = 2 -> "during the second year" -> reduce number to year 1; + 1 to be between week 1 and 52
  set dist_inf 0
  set week_inf_max 0
  set list_repro [ 0.06    ;; January
                   0.16    ;; February
                   0.39    ;; March
                   0.73    ;; April
                   0.8     ;; May
                   0.88    ;; June
                   0.94    ;; July
                   0.97    ;; August
                   1.0     ;; September
                   0.0     ;; October
                   0.0     ;; November
                   0.0 ]   ;; December

  ;; INPUT DATA
  file-open "RadiusAll.txt"
  set list_n_max []
  while [ not file-at-end? ] [ set list_n_max lput file-read list_n_max ]
  file-close

  ;; OUTPUT
  set list_inf_struc []
  set list_trans_struc []
  set list_leth_struc []
  set list_mue_curr []
  set list_mue_all []

end


;; LOAD LANDSCAPE _____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
  ;; observer procedure
  ;; Load a landscape as .txt file.
  ;; Create a "standardised outbreak area" with mean quality.
  ;; The darker the color is, the higher is the breeding capacity of the home range.

to Init_Landscape

  ;; file import
  ifelse (file = "" OR file = "hom")
  [
    ask patches [ set quality mean_quality ]
  ][
    ifelse seed_setup = "BS"
    [
      file-open (word "nlms/" file "_" (behaviorspace-run-number mod rep + 1) ".txt")
    ][
      file-open (word "nlms/" file ".txt")
    ]

    let x 0
    let y 0
    let patch_id 0
    ask patches [ set quality -999]

    while [ not file-at-end? ]
    [
      while [ y < max-pycor + 1 ]
      [
        set x 0
        while [ x < max-pxcor + 1 ]
        [
          let next-value file-read  ;; iterate over file entries seperated by white-space
          ask patch x y
          [
            set quality next-value
            if patch_id = 0
            [
              set id patch_id
              set patch_id patch_id + 1
            ]
            ifelse any? neighbors with [quality = [quality ] of myself AND id > 0]
            [
              let min_id min [id] of neighbors with [quality = [quality] of myself AND id > 0]
              set id min_id
              ask neighbors with [quality = [quality] of myself AND id > min_id] [ set id min_id ]
            ]
            [
              set id patch_id set patch_id patch_id + 1
          ] ]
          set x x + 1
        ]
        set y y + 1
    ] ]
    file-close
  ]

  ; define start patches (borrowing variables habitat and matrix)
  ;set habitat patches with [(pycor = max-pycor OR pycor = max-pycor - 1) AND (pxcor = max-pxcor / 2 OR pxcor = max-pxcor / 2 - 1 OR pxcor = max-pxcor / 2 + 1)]
  set habitat patches with [ pxcor = (max-pxcor / 2) AND pycor = max-pycor ]
  set matrix patches with [not member? self habitat AND quality > 0]

  ask habitat [ set quality mean_quality ]

  let quality_sum sum [quality] of patches
  let quality_sum_max count patches * mean_quality
  let max_quality 2 * mean_quality

  while [ quality_sum > quality_sum_max ]
  [
    let id_c [id] of one-of matrix with [quality >= 0.01 AND quality < max_quality]
    ask matrix with [id = id_c]
    [
      set quality quality - 0.01
      set quality_sum quality_sum - 0.01
  ] ]

  while [ quality_sum < quality_sum_max ]
  [
    let id_c [id] of one-of matrix with [quality > 0 AND quality <= (max_quality - 0.01)]
    ask matrix with [id = id_c]
    [
      set quality quality + 0.01
      set quality_sum quality_sum + 0.01
  ] ]

  ;; set cluster id
  let x 0
  let y 0
  let patch_id 0

  while [ y < max-pycor + 1 ]
    [
      set x 0
      while [ x < max-pxcor + 1 ]
      [
        ask patch x y
        [
          if patch_id = 0
          [
            set id patch_id
            set patch_id patch_id + 1
          ]
          ifelse any? neighbors with [precision quality 2 = [precision quality 2] of myself AND id > 0]
          [
            let min_id min [id] of neighbors with [precision quality 2 = [precision quality 2] of myself AND id > 0]
            set id min_id
            ask neighbors with [precision quality 2 = [precision quality 2] of myself AND id > min_id] [ set id min_id ]
          ]
          [
            set id patch_id set patch_id patch_id + 1
        ] ]
        set x x + 1
      ]
      set y y + 1
  ]

  set habitat patches with [quality > 0]
  set matrix patches with [quality <= 0]
  ask matrix [ set quality 0 ]
  ask patches [ set quality_int ceiling quality]  ;; e.g. 5 is the range from more than 4 to 5
  let qual_mean_temp mean [quality] of patches

  ask habitat
  [
    set capacity (quality * (20 / qual_mean_temp))   ; quality * (mean number of 20 boars per cell / mean quality) = 4.4445 boars per quality unit
    set count_repro (count turtles-here with [is_female = 1 AND is_mother = 0 AND (age_group = "adult" OR age_group = "Subdult")])
    set disp_cells habitat in-radius dist_disp
    set qual_around sum [quality] of neighbors
    set is_habitat 1
  ]
  ask matrix [ set capacity 0 ]

  ask patch (max-pxcor / 2) max-pycor [ set quality_mean_6 mean [quality] of patches in-radius 6 ]

  set patches_0 count matrix
  set patches_1 count habitat with [quality_int = 1]
  set patches_2 count habitat with [quality_int = 2]
  set patches_3 count habitat with [quality_int = 3]
  set patches_4 count habitat with [quality_int = 4]
  set patches_5 count habitat with [quality_int = 5]
  set patches_6 count habitat with [quality_int = 6]
  set patches_7 count habitat with [quality_int = 7]
  set patches_8 count habitat with [quality_int = 8]
  set patches_9 count habitat with [quality_int = 9]

end


;; CREATE START POPULATION ____________________________________________________________________________________________________________________________________________________________________________________________________________________________________
  ;; observer procedure
  ;; One wild boar and its herd is allocated to each habitat patch, group size is initialised as three times breeding capacity (quality) of each patch.
  ;; Initial age distributions are taken from the results of a 100 years model run conducted by Kramer-Schadt et al. 2009.

to Init_Population

  crt (herd_prop / 100 * count habitat)
  [
     set age 0
     set dem_stat "dsStart"
     set epi_stat "esSusc"
     set is_shedding 0
     set is_breeder 0
     set is_mother 0
  ]

  ask turtles [ ht ]

  let temp_herd_id 1

  while [ any? turtles with [dem_stat = "dsStart"] ]
  [
    ask turtles with [herd_id = 0]
    [
      move-to one-of habitat with [is_occupied = 0]
      set temp_herd_id temp_herd_id + 1
      set herd_id temp_herd_id
      set dem_stat "dsResident"
      let offspring (release_fct * quality - 1)
      hatch offspring
      if random-float 1 < (offspring - floor(offspring)) [ hatch 1 ]
      set is_occupied 1 ;patches-own
  ] ]

  while [ any? turtles with [age = 0] ]
  [
    ask turtles with [age = 0]
    [
      ifelse random-float 1 <= fem_prob [ set is_female 1 ] [ set is_female 0 ]
      let rand random-float 1
      ;; generate variation in age of start population
      let Var round random-normal 0 (age_blur / 2)
      if Var <= age_blur AND Var >= (- age_blur)
      [
        ;; intial age distribution based on Kramer-Schadt et al. 2009
        if rand <= 0.38                     [set age  52 + Var]
        if rand > 0.38 AND rand <= 0.62 [set age 104 + Var]
        if rand > 0.62 AND rand <= 0.77 [set age 156 + Var]
        if rand > 0.77 AND rand <= 0.86 [set age 208 + Var]
        if rand > 0.86 AND rand <= 0.92 [set age 260 + Var]
        if rand > 0.92 AND rand <= 0.95 [set age 312 + Var]
        if rand > 0.95 AND rand <= 0.97 [set age 364 + Var]
        if rand > 0.97 AND rand <= 0.98 [set age 416 + Var]
        if rand > 0.98 AND rand <= 0.99 [set age 468 + Var]
        if rand > 0.99 AND rand <= 1    [set age 520 + Var]
      ]

      ;; age classes males
      if age <= 21 AND is_female = 0 [ set age_group "piglet" ]
      if age > 21 AND age <= 104 AND is_female = 0 [ set age_group "subadult" ]
      if age > 104 AND is_female = 0
      [
        set age_group "adult"
        set dem_stat "dsRoaming"
      ]

      ;; age classes females
      if age <= 34 AND is_female = 1 [ set age_group "piglet" ]
      if age > 34 AND age <= 52 AND is_female = 1 [ set age_group "subadult" ]
      if age > 52 AND is_female = 1 [ set age_group "adult" ]
  ] ]

;; OUTPUT
  set list_pop_struc []
  ask turtles [ set list_pop_struc fput age list_pop_struc ]
  set list_pop_struc map [ ?1 -> ?1 / 52 ] list_pop_struc

end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                                                                                                                                                                                                                           ;;
;;  GO PROCEDURES                                                                                                                                                                                                                                            ;;
;;                                                                                                                                                                                                                                                           ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;  DISEASE PROCEDURES                                                                                                                                                                                                                                       ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;; DISEASE RELEASE ____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
  ;; observer procedure
  ;; One wild boar group (cell) in the middle of the top row gets lethally infected in a random week ("week_release") of a given year (based on input "year_release").

to Disease_Release

  set tick_leth_max ticks + (mue_max * mue)

  ask patch (max-pxcor / 2) max-pycor
  [
    set inf_x pxcor
    set inf_y pycor
    set count_init count turtles-here
    set count_init_roaming count turtles-here with [ dem_stat = "dsRoaming" ]
    ask turtles-here [ Disease_Infection ]
  ]

end


;; DISEASE TRANSMISSION _______________________________________________________________________________________________________________________________________________________________________________________________________________________________________
  ;; observer procedure
  ;; First, infection pressure is calculated for all habitat cells based on infected indivudals in the herd and the eight herds surrounding the home range cell.
  ;; Afterwards, the model iterates over all individuals and stochastically sets susceptible indivudals to infected if a drawn random number (0-1) is smaller than the infection pressure of the cell.
  ;; The transmission parameter was reversly fitted to the estimated disease spread velocity of approx. 8 km per quarter (Rossi et al. 2010). The transmission probabilities were assigned constant as 0.0208 within a herd and 0.00208 between herds.
  ;; Based on the case fatality rate (case_fatality), on infection the host is stochastically assigned either as lethally infected (esLeth) or transiently infected (esTrans).
  ;; case_fatality for lethally infected individuals applies unchanged for yearlings (subadults), is decreased for adults (case_fatality ^ 2) and increased for piglets (sqrt case_fatality) to represent age-dependent disease outcomes (Dahle & Liess 1992).
  ;; Lethally infected hosts receive their individual infectious period / survival times (in weeks) which are drawn from an exponential distribution with mean mue (resulting in tick_leth).

to Disease_Transmission

  set new_leth 0
  set new_trans 0
  set trans_g 0
  set tick_leth_max ticks + (mue_max * mue)

  ;; update infection status of cells
  ask habitat [ set is_infectious 0 ]

  ask habitat
  [
    ifelse roaming = "OFF"
    [ set count_inf count turtles-here with [is_shedding = 1] ]   ; all individuals are equal -> count all infectious individuals per cell
    [ set count_inf count turtles-here with [(dem_stat != "dsRoaming") AND is_shedding = 1] ]   ; age- and sex-dependent differences -> count only infectious group mates (not adult males)
  ]

  ask habitat
  [
    ifelse roaming = "OFF"
    [ set inf_press (1 - (1 - b_within) ^ count_inf * (1 - b_between) ^ sum [count_inf] of neighbors) ]   ; all individuals are equal -> focal cell and neighboring cells
    [ set inf_press (1 - (1 - b_within) ^ count_inf) ]   ; age- and sex-dependent differences -> only the focal cell
  ]

  ask turtles with [dem_stat != "dsRoaming" AND epi_stat = "esSusc" ]
  [
    if random-float 1 < [inf_press] of patch-here
    [
      Disease_Infection
      set trans_g trans_g + 1
  ] ]

end


;; DISEASE STATE TRANSITION ___________________________________________________________________________________________________________________________________________________________________________________________________________________________________
  ;; observer procedure
  ;; Transient shedders (esTrans) are converted to immune (esImm) after t_trans.
  ;; Individuals protected by maternal antibodies (esImmMat) turn suscpetible (esSusc) after reaching the maximum duration of immunity of 12 weeks.

to Disease_Transition

  ;; transient infecteds
  ask turtles with [epi_stat = "esTrans" AND tick_inf + t_trans <= ticks]
  [
    set epi_stat "esImm"
    set is_shedding 0
  ]
  ;; maternally protected pigs
  ask turtles with [epi_stat = "esImmMat" AND age >= t_anti]
  [
    set epi_stat "esSusc"
    set is_matAB 0
  ]

  ask turtles with [is_matAB = 1 AND age > t_anti - 4]   ;; 1 month with varying immunity for piglets
  [
    if random-float 1 < 0.5
    [
      set epi_stat "esSusc"
      set is_matAB 0
  ] ]

  ;; OUTPUT
  set list_inf_struc []
  ask turtles with [is_shedding = 1]  [ set list_inf_struc fput age list_inf_struc ]
  set list_inf_struc map [ ?1 -> ?1 / 52 ] list_inf_struc

  set list_trans_struc []
  ask turtles with [epi_stat = "esTrans"]  [ set list_trans_struc fput age list_trans_struc ]
  set list_trans_struc map [ ?1 -> ?1 / 52 ] list_trans_struc

  set list_leth_struc []
  ask turtles with [epi_stat = "esLeth"]  [ set list_leth_struc fput age list_leth_struc ]
  set list_leth_struc map [ ?1 -> ?1 / 52 ] list_leth_struc

end


;; INFECTIOUS PROCESS _________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
  ;; turtle procedure
  ;; Individual infection based on case fatality

to Disease_Infection

  if epi_stat = "esSusc" AND age_group = "piglet"   [ set ind_case_fatality sqrt case_fatality ]
  if epi_stat = "esSusc" AND age_group = "subadult" [ set ind_case_fatality case_fatality ]
  if epi_stat = "esSusc" AND age_group = "adult"    [ set ind_case_fatality case_fatality ^ 2 ]

  ifelse random-float 1 <= ind_case_fatality
  [
    set epi_stat "esLeth"
    set is_shedding 1
    set tick_inf ticks + 1
    set tick_leth 999   ;; arbitrary high value
    while [ tick_leth > tick_leth_max ] [ set tick_leth ticks + 1 + floor (- (mue - 0.5) * ln random-float 1) ]
    set list_mue_all fput (tick_leth - ticks) list_mue_all  ;; list for lethally infected individuals only
    set new_leth new_leth + 1
    set is_infectious 1  ;patches-own
    set is_infected 1  ;patches-own
    if distancexy inf_x inf_y > dist_inf [ set dist_inf distancexy inf_x inf_y ]  ;patches-own
    set inf_id id
  ][ ; else
    set epi_stat "esTrans"
    set is_shedding 1
    set tick_inf ticks + 1
    set new_trans new_trans + 1
    set is_infectious 1  ;patches-own
    set is_infected 1  ;patches-own
    if distancexy inf_x inf_y > dist_inf [ set dist_inf distancexy inf_x inf_y ]  ;patches-own
    set inf_id id
  ]

end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;  HOST PROCEDURES                                                                                                                                                                                                                                          ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;; REPRODUCTION _______________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
  ;; observer procedure
  ;; Adult and subadult females reproduce only once a year depending on their age class.
  ;; In the first week of the year females are checked wether they are able breed (age + no piglets) and if there is at least one reproducible male within their range (> 2 years).
  ;; Litter size is drawn from a truncated normal distribution and reduced to a constant fraction for infected individuals.
  ;; Depending on the disease state of breeding indiviuals its piglets disease states are adjusted.

to Host_Reproduction

  ask habitat [ set count_repro count turtles-here with [is_female = 1 AND is_mother = 0 AND age > 34] ]

;  ;; at the beginning of each year check each cell for reproductive females (no age-based hierarchy)
;  if week = 1
;  [
;    ask habitat with [count_repro > 0]
;    [
;      let repros turtles-here with [is_female = 1 AND is_mother = 0 AND age > 34]
;      ifelse count_repro > quality
;      [
;        ask n-of quality repros
;        [
;          set is_breeder 1
;          set breed_prob random-float 1
;      ] ]
;      [  ;; else
;        ask repros
;        [
;          set is_breeder 1
;          set breed_prob random-float 1
;  ] ] ] ]

  ;; at the beginning of each year check... (age-based hierarchy)
  if week = 1
  [
    ;; ... for habitat patches with at least one male to reproduce ...
    ask habitat
    [
      ;; ... and for all females which can breed
      let fems count_repro
      ;; if number of potential breeders is higher than breeding capacity
      let qual_rand quality_int
      if random-float 1 < (quality - floor(quality)) [ set qual_rand floor quality ]
      if (fems > qual_rand)
      [
        ;; … sort breeders descending by age …
        let fems_sort sort-on [(- age)] turtles-here with [is_female = 1 AND is_mother = 0 AND age > 34]
        set fems qual_rand
        if length fems_sort > 0
        [
          while [ fems > 0 ]
          [
            ;; … and let num-fem of them (difference between capacity and number of breeders) reproduce
            ask first fems_sort
            [
              set is_breeder 1
              set breed_prob random-float 1
            ]
            set fems_sort remove-item 0 fems_sort
            set fems (fems - 1)
      ] ] ]

      ;; if number of breeders is equal or below than breeding capacity all breeders reproduce
      if (fems > 0 AND fems <= qual_rand)
      [
        ask turtles-here with [is_female = 1 AND age_group = "adult" AND is_mother = 0]
        [
          set is_breeder 1
          set breed_prob random-float 1
  ] ] ] ]

  ;; cumulative reproduction probabilities for each month (January to September)
  ifelse week <= 4 [ set repro_prob item 0 list_repro  Host_ProbRepro ][
    ifelse week <= 8 [ set repro_prob item 1 list_repro  Host_ProbRepro ][
      ifelse week <= 13 [ set repro_prob item 2 list_repro  Host_ProbRepro ][
        ifelse week <= 17 [ set repro_prob item 3 list_repro  Host_ProbRepro ][
          ifelse week <= 21 [ set repro_prob item 4 list_repro  Host_ProbRepro ][
            ifelse week <= 26 [ set repro_prob item 5 list_repro  Host_ProbRepro ][
              ifelse week <= 30 [ set repro_prob item 6 list_repro  Host_ProbRepro ][
                ifelse week <= 34 [ set repro_prob item 7 list_repro  Host_ProbRepro ][
                  ifelse week <= 39 [ set repro_prob item 8 list_repro  Host_ProbRepro ][
                    ifelse week <= 43 [ set repro_prob item 9 list_repro  Host_ProbRepro ][
                      ifelse week <= 47 [ set repro_prob item 10 list_repro  Host_ProbRepro ][
                        if week <= 52 [ set repro_prob item 11 list_repro  Host_ProbRepro
  ] ] ] ] ] ] ] ] ] ] ] ]

end

to Host_ProbRepro
  ;; observer procedure

    ask turtles with [is_breeder = 1]
    [
      ;; stochastic effect in probabilities of reproduction
      if breed_prob < repro_prob
      [
        let rand random-float 1

        ;; stochastic effect in number of piglets
        if (rand > 0.01305859 AND rand <= 0.082208805)
        [
          ;; if the breeder is infected number of piglets is reduced
          ifelse is_shedding = 1
          [
            hatch floor (1 * fert_red)   ;; round to integer
            [
              set dem_stat "dsOffspring"
              set family who
          ] ]
          [
          ;; if the breeder is not infected
            hatch 1
            [
              set dem_stat "dsOffspring"
              set family who
        ] ] ]

        if (rand > 0.082208805 AND rand <= 0.24510931)
        [
          ;; if the breeder is infected number of piglets is reduced
          ifelse is_shedding = 1
          [
            hatch floor (2 * fert_red)
            [
              set dem_stat "dsOffspring"
              set family who
            ]
          ][
            ;; if the breeder is not infected
            hatch 2
            [
              set dem_stat "dsOffspring"
              set family who
        ] ] ]

        if (rand > 0.24510931 AND rand <= 0.495050085)
        [
          ;; if the breeder is infected number of piglets is reduced
          ifelse is_shedding = 1
          [
            hatch floor (3 * fert_red)
            [
              set dem_stat "dsOffspring"
              set family who
            ]
          ][
            ;; if the breeder is not infected
            hatch 3
            [
              set dem_stat "dsOffspring"
              set family who
        ] ] ]

        if (rand > 0.495050085 AND rand <= 0.744990859)
        [
          ;; if the breeder is infected number of piglets is reduced
          ifelse is_shedding = 1
          [
            hatch floor (4 * fert_red)
            [
              set dem_stat "dsOffspring"
              set family who
          ] ]
          [
            ;; if the breeder is not infected
            hatch 4
            [
              set dem_stat "dsOffspring"
              set family who
        ] ] ]

        if (rand > 0.744990859 AND rand <= 0.907891364)
        [
          ;; if the breeder is infected number of piglets is reduced
          ifelse is_shedding = 1
          [
            hatch floor (5 * fert_red)
            [
              set dem_stat "dsOffspring"
              set family who
            ]
          ][
            ;; if the breeder is not infected
            hatch 5
            [
              set dem_stat "dsOffspring"
              set family who
        ] ] ]

        if (rand > 0.907891364 AND rand <= 0.977041579)
        [
          ;; if the breeder is infected number of piglets is reduced
          ifelse is_shedding = 1
          [
            hatch floor (6 * fert_red)
            [
              set dem_stat "dsOffspring"
              set family who
            ]
          ][
            ;; if the breeder is not infected
            hatch 6
            [
              set dem_stat "dsOffspring"
              set family who
        ] ] ]

        if (rand > 0.977041579 AND rand <= 0.996144138)
        [
          ;; if the breeder is infected number of piglets is reduced
          ifelse is_shedding = 1
          [
            hatch floor (7 * fert_red)
            [
              set dem_stat "dsOffspring"
              set family who
          ] ]
          [
            ;; if the breeder is not infected
            hatch 7
            [
              set dem_stat "dsOffspring"
              set family who
        ] ] ]

        if (rand > 0.996144138 AND rand <= 0.999575749)
        [
          ;; if the breeder is infected number of piglets is reduced
          ifelse is_shedding = 1
          [
            hatch floor (8 * fert_red)
            [
              set dem_stat "dsOffspring"
              set family who
          ] ]
          [
            ;; if the breeder is not infected
            hatch 8
            [
              set dem_stat "dsOffspring"
              set family who
        ] ] ]

        if (rand > 0.999575749 AND rand <= 0.99997575)
        [
          ;; if the breeder is infected number of piglets is reduced
          ifelse is_shedding = 1
          [
            hatch floor (9 * fert_red)
            [
              set dem_stat "dsOffspring"
              set family who
          ] ]
          [
            ;; if the breeder is not infected
            hatch 9
            [
              set dem_stat "dsOffspring"
              set family who
        ] ] ]

        if (rand > 0.99997575 AND rand <= 1.0)
        [
          ;; if the breeder is infected number of piglets is reduced
          ifelse is_shedding = 1
          [
            hatch floor (10 * fert_red)
            [
              set dem_stat "dsOffspring"
              set family who
          ] ]
          [
            ;; if the breeder is not infected
            hatch 10
            [
              set dem_stat "dsOffspring"
              set family who
        ] ] ]

        set tick_birth ticks
        set is_breeder 0
        set is_mother 1
  ] ]

  ask turtles with [dem_stat = "dsOffspring"]
  [
    ifelse random-float 1 <= fem_prob [ set is_female 1 ] [ set is_female 0 ]
    set age 0
    set age_group "piglet"
    set is_breeder 0

    ;; Susceptibles
    if [epi_stat] of turtle family = "esSusc"
    [
      set dem_stat "dsResident"
      set epi_stat "esSusc"
      set is_shedding 0
    ]

    ;; Immunes
    if [epi_stat] of turtle family = "esImm"
    [
      set dem_stat "dsResident"
      set epi_stat "esImmMat"
      set is_shedding 0
      set is_matAB 1
    ]

    ;; Infecteds
    if [is_shedding] of turtle family = 1
    [
      ifelse random-float 1 < fetal_inf
      [
        set dem_stat "dsResident"
        set epi_stat "esLeth"
        set is_shedding 1
        set tick_inf ticks + 1
        set tick_leth 999   ;; arbitrary high value
        while [ tick_leth > tick_leth_max ] [ set tick_leth ticks + 1 + floor (- (mue - 0.5) * ln random-float 1) ]
        set new_leth new_leth + 1
        set list_mue_all fput (tick_leth - ticks) list_mue_all
      ][
        set dem_stat "dsResident"
        set epi_stat "esSusc"
        set is_shedding 0
  ] ] ]

end


;; DEATH ______________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
  ;; observer procedure
  ;; Iterating over the entire population, each individual either stochastically dies due to age class dependent mortality rates by drawing a random number,
  ;;    due to reaching a certain maximum age or due to a lethal infection after a certain ifnection time span (tick_leth).
  ;; Stochastic baseline mortality is age-dependent and are drawn from a Gaussian distribution which was adjusted to annual survival estimates found in the literature and is cutted symmetrically around the mean.
  ;; The stochastic effect resembles 'good' versus 'bad' years for wild boars (i.e. environmental noise).
  ;; Each time step the adjusted age-dependent mortality is applied to all individual.


to Host_Death

  ;; environmental noise (turned off)
  ;if week = 1
  ;[
  ;  let Noise 999
  ;  while [ Noise > 1 OR Noise < -1 ]
  ;  [
  ;    set Noise random-normal 0 (1 / 2.576)   ; 99% between -1 and 1
  ;    ;; survival probability for subadults and adults
  ;    set survival_prob 0.65 + Noise * (0.65 - 0.4)   ; mean + Noise * (mean - min) to scale for age group
  ;    ;; survival probability for piglets
  ;    set survival_prob_pig 0.5 + Noise * (0.5 - 0.1)
  ;  ]
  ;]

  ;; survival probability for subadults and adults
  set survival_prob 0.65
  ;; survival probability for piglets
  set survival_prob_pig 0.5

  ask turtles
  [
    ifelse age > 34
    [
      if random-float 1 < (1 - (survival_prob ^ (1 / 52))) [ die ]
    ][  ;; mean survival probability for yearlings (Gaillard et al. 1987) = 0.65 = mean survival probability for adults (Focardi et al. 1996)
      if random-float 1 < (1 - (survival_prob_pig ^ (1 / 52))) [ die ]
  ] ]

  ;; mortality due to a lethal infection
  ask turtles with [epi_stat = "esLeth"]
  [
    if tick_leth = ticks [ die ]
  ]

  ;; mortality due to reaching maximum age (more than 10 years)
  ask turtles with [age >= max_age] [ die ]

  ask habitat with [not any? turtles-here AND is_occupied = 1] [ set is_occupied 0 ]

  ;;OUTPUT
  set list_pop_struc []
  ask turtles [ set list_pop_struc fput age list_pop_struc ]
  set list_pop_struc map [ ?1 -> ?1 / 52 ] list_pop_struc

end


;; AGE ________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
  ;; observer procedure
  ;; For each individual age is incremented one week and disease state transitions are performed.
  ;; Females who gave birth more than 33 weeks before can reproduce again.

to Host_Ageing

  ask turtles
  [
    set age age + 1
    ifelse age = 35
    [
      set age_group "subadult"
      if is_female = 0 [ set dem_stat "dsDispM_subadult" ]
    ][
      ifelse age = 53 AND is_female = 1
      [ set age_group "adult" ]
      [
        if age = 105 AND is_female = 0
        [
          set age_group "adult"
          if qual_around > 0 [ set dem_stat "dsRoaming" ]  ;patches-own
    ] ] ]
    if is_mother = 1 AND (tick_birth + 34) = ticks [ set is_mother 0 ]
  ]

end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;  HOST MOVEMENT PROCEDURES                                                                                                                                                                                                                                 ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;; HERD SPLITTING _____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
  ;; observer procedure
  ;; During summer, subadult females either remain with their mothers or establish new territories nearby.
  ;; Herd splitting is performed in one week per year only where female subadults (8 months = 34 weeks < age <= 2 years = 104 weeks) without offspring are able to move between herds.
  ;; All herds to split are extracted matching the conditions of containing a number of females exceeding the cells carrying capacity and containing at least 2 dispersers.
  ;; For each splittable herd an habitat cell without female groups and a carrying capacity above 0 within the dispersal distance of 6 km (3 cells) is selected randomly.
  ;; All dispersers from the source cell find a new herd on that cell.

to Movement_DispFemales

  let overcap habitat with [count turtles-here with [is_female = 1 AND age > 34] - quality >= 2 AND qual_around > 0]

  ask overcap
  [
    ;; count females above breeding capacity (quality)
    let num-overfem ceiling (count turtles-here with [is_female = 1 AND age > 34] - quality)
    ;; ... assign those females a disperser-status...
    let num-disp count turtles-here with [age_group = "subadult" AND is_female = 1 AND is_breeder = 0]
    ;; ... which are over breeding capacity ...
    if num-disp > num-overfem [ set num-disp num-overfem ]
    ;; ... and if there are at least 2 dispersers ...
    if num-disp >= 2
    [ ask n-of num-disp turtles-here with [age_group = "subadult" AND is_female = 1 AND is_breeder = 0] [ set dem_stat "dsDispF" ] ]
    set count_disp_fem 0
  ]

  ask overcap with [any? turtles-here with [dem_stat = "dsDispF"]]
  [
    let list_disp_f []
    ask turtles-here with [dem_stat = "dsDispF"] [ set list_disp_f fput self list_disp_f ]
    set count_disp_fem length list_disp_f   ; number of dispersing females of this cell
    let temp_disp_cells disp_cells with [not any? turtles-here with [is_female = 1]]
    ifelse count temp_disp_cells > 0
    [
      let newhabitat one-of temp_disp_cells
      while [ not empty? list_disp_f ]
      [
        ask first list_disp_f
        [
          move-to newhabitat
          set dem_stat "esResident"
          if is_shedding = 1
          [
            set is_infectious 1 ;patches-own
            set is_infected 1 ;patches-own
            if distancexy inf_x inf_y > dist_inf [ set dist_inf distancexy inf_x inf_y ] ;patches-own
          ]
          set list_disp_f remove-item 0 list_disp_f
      ] ]
      ask newhabitat [ set is_occupied 1 ]
      set list_disp_f []
    ][
      ;; else: version used in Fernandez et al. 2006 -> dispersers stay in their maternal cell if no empty habitat cell is available
      ask turtles-here with [dem_stat = "dsDispF"] [ set dem_stat "esResident" ]
  ] ]

end


;;; MALE DISPERSAL _____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
  ;; observer procedure
  ;; Males perform their natal dispersal in an age of 8 months (35 weeks) upwards.
  ;; All dispersing males form a group of at least 3 subadults and move to a cells below capacity with the highest habitat quality in dispersal distance.

to Movement_DispMales

  let natal_males 0

  ask habitat with [any? turtles-here with [dem_stat = "dsDispM_subadult"] AND qual_around > 0]
  [
    set natal_males turtles-here with [dem_stat = "dsDispM_subadult"]
    if count natal_males > 2
    [
      let new_home disp_cells with [count turtles-here < capacity]
      if count new_home > 0
      [
        set new_home max-one-of new_home [quality]
        ask natal_males
        [
          move-to new_home
          set dem_stat "dsResident"
  ] ] ] ]

end


;; RANGING BEHAVIOUR __________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
  ;; observer procedure
  ;; Long-distance movement of males older than 2 years (104 weeks) following the 'ranging strategy'.
  ;; The ranging strategy, usually employed by males at night, consists of ranging outside the activity zone (Morelle et al. 2015).
  ;;
  ;; Ranging individuals move until the individuals’ weekly roaming distance (dist_move) is reached or there was no better decision to make.
  ;; Transmission probabilities on the move (b_move) are re-scaled to the individual movement distance (ind_b_move = b_move / dist_move)
  ;;
  ;; male dispersal distances (Morelle et al. 2015):
  ;;   - daily range covered = 45% of annual range (1.3 to 2.4 km^2)
  ;;   - mean daily distance travelled = 7.2 to 11.4 km
  ;;   - daily moved distance = 3 to 4 km with a maximum of 12 km  ->  per week: 21-28 km (10.5-14 cells), max. 84 km (42 cells)

to Movement_Roaming

  set trans_w 0
  set trans_b 0

  set contacts_1 0
  set contacts_2 0
  set contacts_3 0
  set contacts_4 0
  set contacts_5 0
  set contacts_6 0
  set contacts_7 0
  set contacts_8 0
  set contacts_9 0

  ;ask turtles with [tag = 1]
  ;[
  ;  set tcolor tcolor - 0.05
  ;  set color tcolor
  ;]

  if roaming != "OFF"
  [
    ask patches
    [
      set count_group count turtles-here with [dem_stat != "dsRoaming"]
      set count_male count turtles-here with [dem_stat = "dsRoaming"]
      set count_repro count turtles-here with [is_female = 1 AND age > 34]
      set scent 0
    ]

    ask turtles with [dem_stat = "dsRoaming"]
    [
      set contacts 0
      let dist 999   ;; arbitrary high value
      while [ dist > 42 ] [ set dist (Report_Weibull 1.3 26) + 1 ]
      set dist_move dist
      set ind_b_move (b_move / dist_move)
      let list_visits []

  ;;-------------------------------------------------------------------------------------------------------------------
  ;; MOVEMENT ALGORITHM 1 (random walk)
      if roaming = "RW"
      [
        while [ dist > 0 ]
        [
          let aim neighbors with [is_habitat = 1]
          if count aim = 0
          [
            set list_visits fput patch-here list_visits
            set scent scent + (dist / dist_move)  ;patches-own
            if is_shedding = 1
            [
              repeat dist
              [
                Disease_Shedding_Move
                Movement_Contacts
            ] ]
            if (count_inf > 0 AND epi_stat = "esSusc") [ Disease_Infection_Move ]
            stop
          ]
          move-to one-of aim
          set dist dist - 1
          set list_visits fput patch-here list_visits
          set scent scent + (1 / dist_move)  ;patches-own
          if is_shedding = 1
          [
            Disease_Shedding_Move
            Movement_Contacts
          ]
          if (count_inf > 0 AND epi_stat = "esSusc") [ Disease_Infection_Move ]
      ] ]

  ;;-------------------------------------------------------------------------------------------------------------------
  ;; MOVEMENT ALGORITHM 2 (correlated random walk)
      if roaming = "CRW"
      [
        while [ dist > 0 ]
        [
          let aim neighbors with [is_habitat = 1]
          if count aim = 0
          [
            set list_visits fput patch-here list_visits
            set scent scent + (dist / dist_move)  ;patches-own
            if is_shedding = 1
            [
              repeat dist
              [
                Disease_Shedding_Move
                Movement_Contacts
            ] ]
            if (count_inf > 0 AND epi_stat = "esSusc") [ Disease_Infection_Move ]
            stop
          ]
          ;; wrapped cauchy distribution
          ;; (equation from Haefner & Crist 1994, also used by Zollner & Lima 1999, Fletcher 2006)
          ;; produces values between 0-180 and 540-720 if q = 0
          rt (2 * (atan (( (1 - q) / (1 + q)) * tan ((random-float 1 - 0.5) * 180)) 1))
          while [ patch-ahead 1 = nobody OR not member? patch-ahead 1 aim ] [ rt 30]
          move-to patch-ahead 1
          set dist dist - 1
          set list_visits fput patch-here list_visits
          set scent scent + (1 / dist_move)  ;patches-own
          if is_shedding = 1
          [
            Disease_Shedding_Move
            Movement_Contacts
          ]
          if (count_inf > 0 AND epi_stat = "esSusc") [ Disease_Infection_Move ]
      ] ]

  ;;-------------------------------------------------------------------------------------------------------------------
  ;; MOVEMENT ALGORITHM 3 (habitat-dependent movement)
      if roaming = "HD"
      [
        while [ dist > 0 ]
        [
          let aim neighbors with [is_habitat = 1]
          if count aim = 0
          [
            set list_visits fput patch-here list_visits
            set scent scent + (dist / dist_move)  ;patches-own
            if is_shedding = 1
            [
              repeat dist
              [
                Disease_Shedding_Move
                Movement_Contacts
            ] ]
            if (count_inf > 0 AND epi_stat = "esSusc") [ Disease_Infection_Move ]
            stop
          ]
          ifelse random-float 1 < q [ set aim one-of aim with-max [quality_int] ]
                                    [ set aim one-of aim ]
          face aim
          move-to aim
          set dist dist - 1
          set list_visits fput patch-here list_visits
          set scent scent + (1 / dist_move)  ;patches-own
          if is_shedding = 1
          [
            Disease_Shedding_Move
            Movement_Contacts
          ]
          if (count_inf > 0 AND epi_stat = "esSusc") [ Disease_Infection_Move ]
      ] ]

  ;;-------------------------------------------------------------------------------------------------------------------
  ;; MOVEMENT ALGORITHM 4 (competition-driven movement)
      if roaming = "CD"
      [
        while [ dist > 0 ]
        [
          let aim neighbors with [is_habitat = 1]
          if count aim = 0
          [
            set list_visits fput patch-here list_visits
            set scent scent + (dist / dist_move)  ;patches-own
            if is_shedding = 1
            [
              repeat dist
              [
                Disease_Shedding_Move
                Movement_Contacts
            ] ]
            if (count_inf > 0 AND epi_stat = "esSusc") [ Disease_Infection_Move ]
            stop
          ]
          set count_male count_male - 1  ;patches-own
          ifelse random-float 1 < q
          [
            ask aim
            [
              set W_M (count_group + count_male) / capacity
              if W_M > 1 [ set W_M 1 ]
            ]
            set aim min-one-of aim [W_M]
            face aim
            move-to aim
          ][
            set aim one-of aim
            face aim
            move-to aim
          ]
          set dist dist - 1
          set list_visits fput patch-here list_visits
          set scent scent + (1 / dist_move)  ;patches-own
          set count_male count_male + 1  ;patches-own
          if is_shedding = 1
          [
            Disease_Shedding_Move
            Movement_Contacts
          ]
          if (count_inf > 0 AND epi_stat = "esSusc") [ Disease_Infection_Move ]
      ] ]

      if roaming = "CDV"
      [
        while [ dist > 0 ]
        [
          if patch-ahead -1 != nobody [ ask patch-ahead -1 [ set last_visited 1 ] ]  ;patches-own
          let aim neighbors with [is_habitat = 1 AND last_visited = 0]
          if patch-ahead -1 != nobody [ ask patch-ahead -1 [ set last_visited 0 ] ]  ;patches-own
          if count aim = 0
          [
            set list_visits fput patch-here list_visits
            set scent scent + (dist / dist_move)  ;patches-own
            if is_shedding = 1
            [
              repeat dist
              [
                Disease_Shedding_Move
                Movement_Contacts
            ] ]
            if (count_inf > 0 AND epi_stat = "esSusc") [ Disease_Infection_Move ]
            stop
          ]
          set count_male count_male - 1  ;patches-own
          ifelse random-float 1 < q
          [
            ask aim
            [
              set W_M (count_group + count_male) / capacity
              if W_M > 1 [ set W_M 1 ]
            ]
            set aim min-one-of aim [W_M]
            face aim
            move-to aim
          ][
            set aim one-of aim
            face aim
            move-to aim
          ]
          set dist dist - 1
          set list_visits fput patch-here list_visits
          set scent scent + (1 / dist_move)  ;patches-own
          set count_male count_male + 1  ;patches-own
          if is_shedding = 1
          [
            Disease_Shedding_Move
            Movement_Contacts
          ]
          if (count_inf > 0 AND epi_stat = "esSusc") [ Disease_Infection_Move ]
      ] ]

  ;;-------------------------------------------------------------------------------------------------------------------
  ;; MOVEMENT ALGORITHM 5 (density-dependent movement)
      if roaming = "DD"
      [
        while [ dist > 0 ]
        [
          let aim neighbors with [is_habitat = 1]
          if count aim = 0
          [
            set list_visits fput patch-here list_visits
            set scent scent + (dist / dist_move)  ;patches-own
            if is_shedding = 1
            [
              repeat dist
              [
                Disease_Shedding_Move
                Movement_Contacts
            ] ]
            if (count_inf > 0 AND epi_stat = "esSusc") [ Disease_Infection_Move ]
            stop
          ]
          set count_male count_male - 1  ;patches-own
          ifelse random-float 1 < q
          [
            ask aim [ set W_M (count_group + count_male) / capacity ]
            move-to max-one-of aim [W_M]
          ][
            move-to one-of aim
          ]
          set dist dist - 1
          set list_visits fput patch-here list_visits
          set scent scent + (1 / dist_move)  ;patches-own
          set count_male count_male + 1  ;patches-own
          if is_shedding = 1
          [
            Disease_Shedding_Move
            Movement_Contacts
          ]
          if (count_inf > 0 AND epi_stat = "esSusc") [ Disease_Infection_Move ]
      ] ]

  ;;-------------------------------------------------------------------------------------------------------------------
  ;; MOVEMENT ALGORITHM 6 (male-dependent movement)
      if roaming = "MD"
      [
        while [ dist > 0 ]
        [
          let aim neighbors with [is_habitat = 1]
          if count aim = 0
          [
            set list_visits fput patch-here list_visits
            set scent scent + (dist / dist_move)  ;patches-own
            if is_shedding = 1
            [
              repeat dist
              [
                Disease_Shedding_Move
                Movement_Contacts
            ] ]
            if (count_inf > 0 AND epi_stat = "esSusc") [ Disease_Infection_Move ]
            stop
          ]
          set count_male count_male - 1  ;patches-own
          ifelse random-float 1 < q
          [
            ask aim [ set W_M count_male / capacity ]
            move-to min-one-of aim [W_M]
          ][
            move-to one-of aim
          ]
          set dist dist - 1
          set list_visits fput patch-here list_visits
          set scent scent + (1 / dist_move)  ;patches-own
          set count_male count_male + 1  ;patches-own
          if is_shedding = 1
          [
            Disease_Shedding_Move
            Movement_Contacts
          ]
          if (count_inf > 0 AND epi_stat = "esSusc") [ Disease_Infection_Move ]
      ] ]

  ;;-------------------------------------------------------------------------------------------------------------------
  ;; MOVEMENT ALGORITHM 7 (female/male-dependent movement)
      if roaming = "FD"
      [
        while [ dist > 0 ]
        [
          let aim neighbors with [is_habitat = 1]
          if count aim = 0
          [
            set list_visits fput patch-here list_visits
            set scent scent + (dist / dist_move)  ;patches-own
            if is_shedding = 1
            [
              repeat dist
              [
                Disease_Shedding_Move
                Movement_Contacts
            ] ]
            if (count_inf > 0 AND epi_stat = "esSusc") [ Disease_Infection_Move ]
            stop
          ]
          set count_male count_male - 1  ;patches-own
          ifelse random-float 1 < q   ; quality of decision
          [
            ifelse any? aim with [count_repro > 0]
            [
              ask aim
              [
                ifelse count_repro > 0
                [
                  set W_M ((count_repro - count_male) / (count_repro + count_male))
                  if W_M > 1 [ set W_M 1 ]
                ][
                  set W_M -1
                ]
              ]
              move-to max-one-of aim [W_M]
            ][
              move-to min-one-of aim [count_male]
            ]
          ][
            move-to one-of aim
          ]
          set dist dist - 1
          set list_visits fput patch-here list_visits
          set scent scent + (1 / dist_move)  ;patches-own
          set count_male count_male + 1  ;patches-own
          if is_shedding = 1
          [
            Disease_Shedding_Move
            Movement_Contacts
          ]
          if (count_inf > 0 AND epi_stat = "esSusc") [ Disease_Infection_Move ]
      ] ]

  ;;-------------------------------------------------------------------------------------------------------------------
  ;; MOVEMENT ALGORITHM 8 (correlated habitat-dependent movement)
      if roaming = "HD-CRW"
      [
        while [ dist > 0 ]
        [
          let aim neighbors with [is_habitat = 1]
          if count aim = 0
          [
            set list_visits fput patch-here list_visits
            set scent scent + (dist / dist_move)  ;patches-own
            if is_shedding = 1
            [
              repeat dist
              [
                Disease_Shedding_Move
                Movement_Contacts
            ] ]
            if (count_inf > 0 AND epi_stat = "esSusc") [ Disease_Infection_Move ]
            stop
          ]
          ifelse random-float 1 < q
          [
            ;; CRW for long-range movements only - fix q because local movement arises from habitat-dependent movement
            rt (2 * (atan (( (1 - 0.9) / (1 + 0.9)) * tan ((random-float 1 - 0.5) * 180)) 1))   ;; produces values between 0-180 and 540-720 if q = 0
            while [ patch-ahead 1 = nobody OR member? patch-ahead 1 matrix ] [ rt random 90 - 45 ]
            move-to patch-ahead 1
          ]
          [
            move-to one-of neighbors with-max [quality_int]
          ]
          set dist dist - 1
          set list_visits fput patch-here list_visits
          set scent scent + (1 / dist_move)  ;patches-own
          if is_shedding = 1
          [
            Disease_Shedding_Move
            Movement_Contacts
          ]
          if (count_inf > 0 AND epi_stat = "esSusc") [ Disease_Infection_Move ]
      ] ]

  ;;-------------------------------------------------------------------------------------------------------------------
  ;; AFTERWARDS

      set list_visits remove-duplicates list_visits
      set visited length list_visits
  ] ]

  let roam_inf turtles with [dem_stat = "dsRoaming" AND is_shedding = 1]
  ifelse count roam_inf > 0
  [
    set contacts_med median [contacts] of roam_inf
    set contacts_avg mean [contacts] of roam_inf
    set contacts_lwr lower-quartile [contacts] of roam_inf
    set contacts_upr upper-quartile [contacts] of roam_inf
    set contacts_max max [contacts] of roam_inf
    if count roam_inf > 2 [ set contacts_var variance [contacts] of roam_inf ]
  ][
    set contacts_med -999
    set contacts_avg -999
    set contacts_lwr -999
    set contacts_upr -999
    set contacts_max -999
    set contacts_var -999
]

end


to Movement_Contacts
;; turtle procedure

  let contacts_c (count turtles-here with [epi_stat = "esSusc"] * (1 / dist_move))
  set contacts contacts + contacts_c
  if quality_int = 1 [ set contacts_1 contacts_1 + contacts_c ]
  if quality_int = 2 [ set contacts_2 contacts_2 + contacts_c ]
  if quality_int = 3 [ set contacts_3 contacts_3 + contacts_c ]
  if quality_int = 4 [ set contacts_4 contacts_4 + contacts_c ]
  if quality_int = 5 [ set contacts_5 contacts_5 + contacts_c ]
  if quality_int = 6 [ set contacts_6 contacts_6 + contacts_c ]
  if quality_int = 7 [ set contacts_7 contacts_7 + contacts_c ]
  if quality_int = 8 [ set contacts_8 contacts_8 + contacts_c ]
  if quality_int = 9 [ set contacts_9 contacts_9 + contacts_c ]

end


to Disease_Shedding_Move

  let spread_move ind_b_move
  ask turtles-here with [epi_stat = "esSusc"]
  [
    if random-float 1 < spread_move
    [
      Disease_Infection
      ifelse [inf_id] of myself = id [ set trans_w trans_w + 1 ]
                                     [ set trans_b trans_b + 1 ]
  ] ]

end


to Disease_Infection_Move

  if random-float 1 < (1 - (1 - ind_b_move) ^ (count_inf))
  [
    Disease_Infection
    set trans_g trans_g + 1
  ]

end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; HELP FUNCTIONS
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;; WEIBULL DISTRIBUTION REPORTER ______________________________________________________________________________________________________________________________________________________________________________________________________________________________
  ;; observer procedure
  ;; Reports a Weibull distribution with given scale (originally "k") and shape (originally "lambda") parameters

to-report Report_Weibull [a-scale a-shape]

  let xWei random-float 1
  let yWei (1 / a-scale)
  let zWei ((a-shape * (-1 * ln(1 - xWei)))^ yWei)
  report zWei

end


;; FRAGMENTATION INDEX ________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
  ;; observer procedure
  ;; Reports a fragmentation index F, originally developed to evaluate the spatial dispersion of a home range (Mitchell and Powell, 2004)

to Report_Fragmentation

  set F_infected 0
  set F_infectious 0

  let Nmax_infected count patches with [is_infected = 1]
  let Nmax_infectious count patches with [is_infected = 1]

  set Nmax_infected item Nmax_infected list_n_max
  set Nmax_infectious item Nmax_infectious list_n_max

  if Nmax_infected > 0 [ set F_infected (1 - ((sum ([count neighbors with [is_infected = 1]] of patches with [is_infected = 1])) / Nmax_infected )) ]
  if Nmax_infectious > 0 [ set F_infectious (1 - ((sum ([count neighbors with [any? turtles-here with [is_shedding = 1]]] of patches with [any? turtles-here with [is_shedding = 1]])) / Nmax_infectious )) ]

end


;; QUARTILES ________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
  ;; observer procedure
  ;; Reports lower and upper quartiles (taken from http://ccl.northwestern.edu/netlogo/models/TrafficBasicAdaptiveIndividuals)

to-report upper-quartile [ xs ]

  let med median xs
  let upper filter [ x -> x > med ] xs
  report ifelse-value (empty? upper) [ med ] [ median upper ]

end

to-report lower-quartile [ xs ]

  let med median xs
  let lower filter [ x -> x < med ] xs
  report ifelse-value (empty? lower) [ med ] [ median lower ]

end


;; UPDATE _____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
  ;; observer procedure
  ;; Update variables (e.g. counts, infection status of patches, output lists) and check stop condition (DONE)

to Update

  ;; update variables
  ask patches
  [
    ifelse any? turtles-here with [is_shedding = 1] [ set is_infectious 1 ]
                                                    [ set is_infectious 0 ]
    set count_male count turtles-here with [dem_stat = "dsRoaming"]
  ]
  ifelse count turtles > 1
  [
    set dens_var variance [ count turtles-here ] of habitat
    set dens_roam_var variance [ scent ] of habitat
  ][
    set dens_var 0
    set dens_roam_var 0
  ]
  ifelse ticks > week_release AND count turtles > 1
  [
    set dens_inf_group_var variance [ count turtles-here with [dem_stat != "dsRoaming" AND is_shedding = 1] ] of habitat
    set dens_inf_roam_var variance [ count turtles-here with [dem_stat = "dsRoaming" AND is_shedding = 1] ] of habitat
  ][
    set dens_inf_group_var 0
    set dens_inf_roam_var 0
  ]

  ;; stop condition
  if ((not any? turtles with [is_shedding = 1]) OR (ticks = run_years * 52)) AND (ticks >= week_release - 1) [ set DONE 1 ]
  if DONE = 1 [ stop ]

  ;; OUTPUT
  set list_mue_curr []
  ask turtles with [epi_stat = "esLeth"] [ set list_mue_curr fput (tick_leth - ticks) list_mue_curr ]
  set week_last ticks + 1
  if week_inf_max = 0 AND (any? patches with [pycor = 0 AND is_infectious = 1]) [ set week_inf_max ticks + 1 ]

end


;; UPDATE VIEW ________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
  ;; observer procedure
  ;; Update GUI by coloring individuals by their healt status (epi_stat) and patches by their quality (grey = non-infectious, orange = infectious, yellow = outbreak start)

to View

  ;ask turtles [ st ]
  ask turtles with [dem_stat = "dsRoaming"] [ st ]
  ask turtles with [epi_stat = "esSusc"] [ set color black ]
  ask turtles with [epi_stat = "esImm" OR epi_stat = "esImmMat"] [ set color 53 ]  ;; green
  ask turtles with [epi_stat = "esLeth"]
  [
    set color 14  ;; red
    set size 0.75
  ]
  ask turtles with [epi_stat = "esTrans"]
  [
    set color 45  ;; yellow
    set size 0.75
  ]
  ask patches [ set pcolor scale-color grey quality 12 0 ]
  ask habitat with [is_infectious = 1] [ set pcolor scale-color orange quality 12 0 ]
  if ticks >= week_release [ ask patch inf_x inf_y [ set pcolor 47 ] ]  ;; pale yellow

end


;; TRACK ______________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
  ;; observer procedure
  ;; Track a subset of infected adult males

to Track

  let num min (list 2 count turtles with [dem_stat = "dsRoaming" AND epi_stat = "esLeth"])
  print num
  ask turtles [ set tag 0 ]
  ask n-of num turtles with [dem_stat = "dsRoaming" AND epi_stat = "esLeth"]
  [
    pd
    set tag 1
    set tcolor 19
    set color tcolor
  ]

end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                                                                                                                                                                                                                           ;;
;; END                                                                                                                                                                                                                                                       ;;
;;                                                                                                                                                                                                                                                           ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
@#$#@#$#@
GRAPHICS-WINDOW
916
11
1291
755
-1
-1
14.7
1
10
1
1
1
0
0
0
1
0
24
0
49
1
1
1
ticks
30.0

PLOT
8
154
353
354
Population
Time [weeks]
Number of Individuals
0.0
10.0
0.0
100.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot count turtles"

PLOT
8
357
353
554
Epidemiological Status (All)
Time [weeks]
Number of Individuals
0.0
10.0
0.0
100.0
true
true
"" "if not any? turtles [ stop ]"
PENS
"Susceptible" 1.0 0 -16777216 true "" "plot count turtles with [epi_stat = \"esSusc\"]"
"Immune" 1.0 0 -13210332 true "" "plot count turtles with [epi_stat = \"esImm\"]"
"Immune Mat" 1.0 0 -10899396 true "" "plot count turtles with [epi_stat = \"esImmMat\"]"
"Lethally Inf" 1.0 0 -8053223 true "" "plot count turtles with [epi_stat = \"esLeth\"]"
"Transient Inf" 1.0 0 -817084 true "" "plot count turtles with [epi_stat = \"esTrans\"]"

PLOT
688
11
911
190
Control
NIL
NIL
0.0
1.0
0.0
1.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" ""

TEXTBOX
866
197
1016
215
1
11
139.9
1

TEXTBOX
692
26
905
185
i
11
139.9
0

BUTTON
695
29
750
62
Setup
Setup\nView
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
752
29
826
62
Run Week
Go\nView
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
828
29
902
62
Run Year
Go\nView
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
811
64
902
97
Run Until End
if DONE = 0 [ Go View ]\n\n;if ticks > WeekRelease [ ask Habitat [ set pcolor scale-color orange inf_count max [inf_count] of patches 0 ]\n;ask Habitat [ set pcolor scale-color cyan NoBoars max [NoBoars] of patches 0 ]\n;ask Habitat with [is_infectious = 1 ] [ set pcolor scale-color orange NoBoars 120 0 ] \n;ask turtles [ ifelse is_shedding = 1 AND DemStat = \"dsDispM_adult\" [ st ] [ ht ] ]\n;if max [scent] of patches > 0 [ ask Habitat [ set pcolor scale-color blue scent max [scent] of patches 0 ] ]
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
695
64
808
97
Run Until Infection
while [ ticks + 1 < (week_release) ] [ Go ]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

PLOT
688
202
911
329
Landscape
NIL
NIL
0.0
1.0
0.0
1.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" ""

TEXTBOX
695
220
906
325
i
11
139.9
0

PLOT
688
340
911
521
Host Population
NIL
NIL
0.0
1.0
0.0
1.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" ""

TEXTBOX
692
357
906
517
i
11
139.9
0

SLIDER
695
361
797
394
herd_prop
herd_prop
0
100
100.0
1
1
NIL
HORIZONTAL

SLIDER
800
361
903
394
release_fct
release_fct
0
10
4.5
0.5
1
NIL
HORIZONTAL

SLIDER
695
397
797
430
longevity
longevity
0
50
11.0
1
1
yrs
HORIZONTAL

SLIDER
800
397
903
430
age_blur
age_blur
0
26
6.0
1
1
weeks
HORIZONTAL

PLOT
689
531
912
755
Disease Dynamics
NIL
NIL
0.0
1.0
0.0
1.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" ""

TEXTBOX
695
549
908
750
i
11
139.9
0

INPUTBOX
698
687
763
747
b_within
0.022
1
0
Number

INPUTBOX
762
687
838
747
b_between
0.0021
1
0
Number

SLIDER
696
551
797
584
year_release
year_release
1
5
2.0
1
1
NIL
HORIZONTAL

SLIDER
696
588
797
621
fetal_inf
fetal_inf
0
1
0.5
0.05
1
NIL
HORIZONTAL

SLIDER
801
551
904
584
fert_red
fert_red
0.025
1
0.625
0.025
1
NIL
HORIZONTAL

INPUTBOX
801
624
851
684
t_trans
1.0
1
0
Number

INPUTBOX
854
624
904
684
t_anti
12.0
1
0
Number

INPUTBOX
695
100
784
160
run_years
12.0
1
0
Number

SLIDER
801
587
904
620
case_fatality
case_fatality
0
1
0.9
0.05
1
NIL
HORIZONTAL

MONITOR
576
107
626
152
Year
ceiling (ticks / 52)
0
1
11

PLOT
357
357
683
554
Age Classes of Infected Ind.
Age Class [years]
Numberof Infected Ind.
0.0
11.0
0.0
100.0
true
true
"" "if not any? turtles with [epi_stat = \"esLeth\" or epi_stat = \"esTrans\"]\n[\n  clear-plot\n  stop\n]"
PENS
"pen-0" 1.0 1 -14737633 true "" ""
"pen-1" 1.0 1 -5298144 true "" ""
"pen-2" 1.0 1 -817084 true "" ""
"pen-3" 1.0 0 -7500403 true "" ""
"oen-4" 1.0 0 -2674135 true "" ""
"Total Inf" 1.0 1 -14737633 true "" "histogram list_inf_struc"
"Lethally Inf" 1.0 1 -5298144 true "" "histogram list_leth_struc"
"Transient Inf" 1.0 1 -817084 true "" "histogram list_trans_struc"

SLIDER
800
433
903
466
dist_disp
dist_disp
0
25
3.0
1
1
NIL
HORIZONTAL

SLIDER
695
433
797
466
fem_prob
fem_prob
0
1
0.5
0.01
1
NIL
HORIZONTAL

TEXTBOX
872
567
1022
585
NIL
11
0.0
1

MONITOR
602
410
676
451
New Transient
new_trans
17
1
10

MONITOR
577
10
681
55
Run Time [min]
runtime
17
1
11

PLOT
357
557
683
756
Survival Time of Lethally Infected Ind.
Duration Infectious Period [weeks]
Number of Infected
0.0
20.0
0.0
10.0
true
false
"" "if not any? turtles with [epi_stat = \"esLeth\"]\n[\n  clear-plot\n  stop\n]"
PENS
"default" 1.0 1 -14737633 true "" "histogram list_mue_curr"

MONITOR
622
578
679
619
max
max(list_mue_all)
17
1
10

PLOT
9
557
354
756
Epidemiological Status (Infected Ind.)
Time [weeks]
Number of Infected Ind.
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"Lethally Inf" 1.0 0 -8053223 true "" "plot count turtles with [epi_stat = \"esLeth\"]"
"Transient Inf" 1.0 0 -817084 true "" "plot count turtles with [epi_stat = \"esTrans\"]"

MONITOR
274
692
351
733
Transient Inf
count turtles with [epi_stat = \"esTrans\"]
17
1
10

MONITOR
274
652
351
693
Lethally Inf
count turtles with [epi_stat = \"esLeth\"]
17
1
10

MONITOR
274
613
351
654
Immunes
count turtles with [epi_stat = \"esImm\" or epi_stat = \"esImmMat\"]
17
1
10

MONITOR
274
574
351
615
Susceptibles
count turtles with [epi_stat = \"esSusc\"]
17
1
10

PLOT
356
154
682
354
Age Classes
Age Class [years]
Number of Individuals
0.0
11.0
0.0
10.0
true
false
"" "if not any? turtles\n[\n  clear-plot\n  stop\n]"
PENS
"default" 1.0 1 -16777216 true "" "histogram list_pop_struc"

MONITOR
558
255
619
296
Males
count turtles with [is_female = 0]
17
1
10

MONITOR
558
215
619
256
Females
count turtles with [is_female = 1]
17
1
10

MONITOR
498
175
559
216
Population
count turtles
17
1
10

MONITOR
615
255
672
296
Adults
count turtles with [age_group = \"adult\"]
17
1
10

MONITOR
615
215
672
256
Subadults
count turtles with [age_group = \"subadult\"]
17
1
10

MONITOR
602
371
676
412
New Lethal
new_leth
17
1
10

MONITOR
631
107
681
152
Week
week
17
1
11

OUTPUT
8
10
572
152
11

MONITOR
577
58
681
103
Ticks per second
ticks / (runtime * 60)
1
1
11

SLIDER
700
285
904
318
mean_quality
mean_quality
1
10
4.5
0.5
1
NIL
HORIZONTAL

BUTTON
787
148
902
181
Default Values
;; CONTROL\nset run_years 12\nset seed_setup \"ON\"\nset seed 1\n\n;; LANDSCAPE\nset file \"hom\"\nset mean_quality 4.5\n\n;; HOST POPULATION\nset herd_prop 100\nset release_fct 4.5\nset longevity 11\nset age_blur 6\nset fem_prob 0.5\nset dist_disp 3\n\n;; DISEASE DYNAMICS\nset year_release 2\nset fert_red 0.625\nset fetal_inf 0.5\nset case_fatality 0.5\nset b_within 0.021\nset b_between 0.0021\nset b_move 0.021\nset mue 4\nset mue_max 10\nset t_trans 1\nset t_anti 12\n\n;; MOVEMENT\nset q 0.7\nset Roaming \"HD\"
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

MONITOR
558
175
619
216
Roamers
count turtles with [dem_stat = \"dsRoaming\"]
17
1
10

SLIDER
801
469
904
502
q
q
0
1.0
0.7
0.05
1
NIL
HORIZONTAL

CHOOSER
696
469
798
514
roaming
roaming
"OFF" "RW" "CRW" "HD" "CD" "CDV" "DD" "MD" "FD" "HD-CRW"
4

CHOOSER
787
100
902
145
seed_setup
seed_setup
"OFF" "ON" "BS"
2

INPUTBOX
836
687
905
747
b_move
0.022
1
0
Number

SLIDER
698
624
797
657
mue
mue
1
9
4.0
1
1
NIL
HORIZONTAL

MONITOR
615
175
672
216
Piglets
count turtles with [age_group = \"piglet\"]
17
1
10

INPUTBOX
734
120
784
180
seed
191.0
1
0
Number

BUTTON
1297
12
1360
45
track
clear-drawing\nTrack
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

INPUTBOX
698
222
904
282
file
patch_l
1
0
String

BUTTON
1297
48
1361
81
export
ask patches [ set pcolor scale-color orange scent 0 5000 ]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
698
651
797
684
mue_max
mue_max
1
15
10.0
1
1
NIL
HORIZONTAL

@#$#@#$#@
## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.0.4
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="Experiment_full" repetitions="200" runMetricsEveryStep="true">
    <setup>Setup</setup>
    <go>Go</go>
    <timeLimit steps="624"/>
    <exitCondition>DONE = 1</exitCondition>
    <metric>seed</metric>
    <metric>count turtles</metric>
    <metric>count turtles with [epi_stat = "esSusc"]</metric>
    <metric>count turtles with [epi_stat = "esTrans"]</metric>
    <metric>count turtles with [epi_stat = "esLeth"]</metric>
    <metric>count turtles with [dem_stat = "dsRoaming"]</metric>
    <metric>count turtles with [dem_stat = "dsRoaming" AND epi_stat = "esSusc"]</metric>
    <metric>count turtles with [dem_stat = "dsRoaming" AND epi_stat = "esTrans"]</metric>
    <metric>count turtles with [dem_stat = "dsRoaming" AND epi_stat = "esLeth"]</metric>
    <metric>new_trans</metric>
    <metric>new_leth</metric>
    <metric>contacts_avg</metric>
    <metric>contacts_med</metric>
    <metric>contacts_lwr</metric>
    <metric>contacts_upr</metric>
    <metric>contacts_max</metric>
    <metric>contacts_var</metric>
    <metric>contacts_1</metric>
    <metric>contacts_2</metric>
    <metric>contacts_3</metric>
    <metric>contacts_4</metric>
    <metric>contacts_5</metric>
    <metric>contacts_6</metric>
    <metric>contacts_7</metric>
    <metric>contacts_8</metric>
    <metric>contacts_9</metric>
    <metric>trans_w</metric>
    <metric>trans_b</metric>
    <metric>trans_g</metric>
    <metric>patches_0</metric>
    <metric>patches_1</metric>
    <metric>patches_2</metric>
    <metric>patches_3</metric>
    <metric>patches_4</metric>
    <metric>patches_5</metric>
    <metric>patches_6</metric>
    <metric>patches_7</metric>
    <metric>patches_8</metric>
    <metric>patches_9</metric>
    <metric>count patches with [is_infected = 1]</metric>
    <metric>count patches with [is_infectious = 1]</metric>
    <metric>max [id] of patches</metric>
    <metric>count_init</metric>
    <metric>count_init_roaming</metric>
    <metric>dens_var</metric>
    <metric>dens_roam_var</metric>
    <metric>dens_inf_group_var</metric>
    <metric>dens_inf_roam_var</metric>
    <metric>quality_mean_6</metric>
    <metric>mean [visited] of turtles with [dem_stat = "dsRoaming"]</metric>
    <metric>F_infectious</metric>
    <metric>F_infected</metric>
    <metric>dist_inf</metric>
    <metric>week_inf_max</metric>
    <metric>week_release</metric>
    <metric>week_last</metric>
    <enumeratedValueSet variable="file">
      <value value="&quot;hom&quot;"/>
      <value value="&quot;rand&quot;"/>
      <value value="&quot;patch_s&quot;"/>
      <value value="&quot;patch_l&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="case_fatality" first="0" step="0.25" last="1"/>
    <steppedValueSet variable="mue" first="1" step="1" last="9"/>
    <enumeratedValueSet variable="mue_max">
      <value value="1"/>
      <value value="3"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="roaming">
      <value value="&quot;RW&quot;"/>
      <value value="&quot;CRW&quot;"/>
      <value value="&quot;HD&quot;"/>
      <value value="&quot;DD&quot;"/>
      <value value="&quot;FD&quot;"/>
      <value value="&quot;MD&quot;"/>
      <value value="&quot;HD-CRW&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="q">
      <value value="0.5"/>
      <value value="0.7"/>
      <value value="0.9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="run_years">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="seed_setup">
      <value value="&quot;BS&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean_quality">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="herd_prop">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="release_fct">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="longevity">
      <value value="11"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="age_blur">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fem_prob">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dist_disp">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="year_release">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fert_red">
      <value value="0.625"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fetal_inf">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t_anti">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t_trans">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="b_within">
      <value value="0.021"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="b_between">
      <value value="0.0021"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="b_move">
      <value value="0.021"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Manuscript_mue4" repetitions="200" runMetricsEveryStep="true">
    <setup>Setup</setup>
    <go>Go</go>
    <timeLimit steps="624"/>
    <exitCondition>DONE = 1</exitCondition>
    <metric>seed</metric>
    <metric>count turtles</metric>
    <metric>count turtles with [epi_stat = "esSusc"]</metric>
    <metric>count turtles with [epi_stat = "esTrans"]</metric>
    <metric>count turtles with [epi_stat = "esLeth"]</metric>
    <metric>count turtles with [dem_stat = "dsRoaming"]</metric>
    <metric>count turtles with [dem_stat = "dsRoaming" AND epi_stat = "esSusc"]</metric>
    <metric>count turtles with [dem_stat = "dsRoaming" AND epi_stat = "esTrans"]</metric>
    <metric>count turtles with [dem_stat = "dsRoaming" AND epi_stat = "esLeth"]</metric>
    <metric>new_trans</metric>
    <metric>new_leth</metric>
    <metric>contacts_avg</metric>
    <metric>contacts_med</metric>
    <metric>contacts_lwr</metric>
    <metric>contacts_upr</metric>
    <metric>contacts_max</metric>
    <metric>contacts_var</metric>
    <metric>contacts_1</metric>
    <metric>contacts_2</metric>
    <metric>contacts_3</metric>
    <metric>contacts_4</metric>
    <metric>contacts_5</metric>
    <metric>contacts_6</metric>
    <metric>contacts_7</metric>
    <metric>contacts_8</metric>
    <metric>contacts_9</metric>
    <metric>trans_w</metric>
    <metric>trans_b</metric>
    <metric>trans_g</metric>
    <metric>patches_0</metric>
    <metric>patches_1</metric>
    <metric>patches_2</metric>
    <metric>patches_3</metric>
    <metric>patches_4</metric>
    <metric>patches_5</metric>
    <metric>patches_6</metric>
    <metric>patches_7</metric>
    <metric>patches_8</metric>
    <metric>patches_9</metric>
    <metric>count patches with [is_infected = 1]</metric>
    <metric>count patches with [is_infectious = 1]</metric>
    <metric>max [id] of patches</metric>
    <metric>count_init</metric>
    <metric>count_init_roaming</metric>
    <metric>dens_var</metric>
    <metric>dens_roam_var</metric>
    <metric>dens_inf_group_var</metric>
    <metric>dens_inf_roam_var</metric>
    <metric>quality_mean_6</metric>
    <metric>mean [visited] of turtles with [dem_stat = "dsRoaming"]</metric>
    <metric>F_infectious</metric>
    <metric>F_infected</metric>
    <metric>dist_inf</metric>
    <metric>week_inf_max</metric>
    <metric>week_release</metric>
    <metric>week_last</metric>
    <enumeratedValueSet variable="case_fatality">
      <value value="0.5"/>
      <value value="0.1"/>
      <value value="0.3"/>
      <value value="0.7"/>
      <value value="0.9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="roaming">
      <value value="&quot;CRW&quot;"/>
      <value value="&quot;HD&quot;"/>
      <value value="&quot;CD&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="file">
      <value value="&quot;hom&quot;"/>
      <value value="&quot;rand&quot;"/>
      <value value="&quot;patch_s&quot;"/>
      <value value="&quot;patch_l&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mue">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mue_max">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="q">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="run_years">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="seed_setup">
      <value value="&quot;BS&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean_quality">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="herd_prop">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="release_fct">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="longevity">
      <value value="11"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="age_blur">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fem_prob">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dist_disp">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="year_release">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fert_red">
      <value value="0.625"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fetal_inf">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t_anti">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t_trans">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="b_within">
      <value value="0.022"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="b_between">
      <value value="0.0021"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="b_move">
      <value value="0.022"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
