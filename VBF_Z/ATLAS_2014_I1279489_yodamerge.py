#! /usr/bin/env python

import yoda, optparse, operator, itertools, sys, math

EFFICIENCIES = ["/ATLAS_2014_I1279489/d01-x01-y02",
                "/ATLAS_2014_I1279489/d01-x02-y02",
                "/ATLAS_2014_I1279489/d02-x01-y02",
                "/ATLAS_2014_I1279489/d02-x02-y02",
                "/ATLAS_2014_I1279489/d01-x01-y03",
                "/ATLAS_2014_I1279489/d01-x02-y03",
                "/ATLAS_2014_I1279489/d02-x01-y03",
                "/ATLAS_2014_I1279489/d02-x02-y03" ]

def MergeEfficiency(path, aos):
  #print "MergeEfficiency: %s"%path
  ## Check that types match, and just output the first one if they don't
  if not all(type(ao) == type(aos[0]) for ao in aos):
    print "WARNING: Several types of analysis object found at path %s: cannot be merged"
    return aos[0]
  ## Check that all of the objects we're about the merge have the same number of points:
  if not all(ao.numPoints == aos[0].numPoints for ao in aos):
    print "WARNING: Not all objects to be merged have the same number of points. Returning first object."
    return aos[0]

  n_points = aos[0].numPoints

  if not all(ao.annotations.has_key("InclusiveSumWeights") for ao in aos):
    print "WARNING: Abandoning merge of path %s because not all inputs have InclusiveSumWeights attribute"%path
    return aos[0]
  
  sumw_tot = sum(float(ao.annotations['InclusiveSumWeights']) for ao in aos)

  point_lists = []
  for i in range(n_points): point_lists.append([])

  for ao in aos:
    counter = 0
    for point in ao.points:
      point_lists[counter].append(point.y*float(ao.annotations['InclusiveSumWeights']))
      counter+=1

  point_counter = 0
  for point in aos[0].points:
    frac = sum(point_lists[point_counter])/sumw_tot
    err = math.sqrt(frac*(1.0-frac)/sumw_tot)
    point.y = frac
    point.yErrs = [err,err]
    point_counter+=1

  aos[0].setAnnotation('InclusiveSumWeights',sumw_tot)
  return aos[0]

def MergeProfile(path, aos):
  print "MergeProfile: %s"%path
  print dir(aos[0]).__iadd__
  ## Check that types match, and just output the first one if they don't
  if not all(type(ao) == type(aos[0]) for ao in aos):
      print "WARNING: Several types of analysis object found at path %s: cannot be merged"
      return aos[0]
  ## Check that all of the objects we're about the merge have the same number of bins:
  if not all(ao.numBins == aos[0].numBins for ao in aos):
    print "WARNING: Not all objects to be merged have the same number of bins. Returning first object."
    return aos[0]

  n_bins = aos[0].numBins
  sumw_tot = sum(float(ao.sumW()) for ao in aos)

  bin_lists = []
  for i in range(n_bins): bin_lists.append([])

  for ao in aos:
    counter = 0
    for bin in ao.bins:
      bin_lists[counter].append(bin.sumWY)
      counter+=1

  print dir(aos[0].bins[0])

  bin_counter = 0
  for bin in aos[0].bins:
    frac = sum(bin_lists[bin_counter])/sumw_tot
    err = math.sqrt(frac*(1.0-frac)/sumw_tot)
    #bin.y = frac
    #bin.yErrs = [err,err]
    bin_counter+=1

  aos[0].setAnnotation('InclusiveSumWeights',sumw_tot)
  return aos[0]

def MergeDistribution(path, aos):
  #print "MergeDistribution: %s"%path
  ## Check that types match, and just output the first one if they don't
  if not all(type(ao) == type(aos[0]) for ao in aos):
      print "WARNING: Several types of analysis object found at path %s: cannot be merged"
      return aos[0]
  ## Check whether normalizations match (possible for histograms via the sumW/integral attrs)
  ## In the absence of better info, we use the norm as a heuristic to change the merging behaviour
  # TODO: What about empty histograms which couldn't be normalized?
  normto = None
  if all(hasattr(ao, "sumW") for ao in aos):
      ## Check that there are some non-empty input histograms
      nonzero_sumws = [ao.sumW() for ao in aos if ao.sumW() != 0]
      sumw_ref = nonzero_sumws[0] if nonzero_sumws else None
      ## Use the first non-empty histogram (if there is one) as a reference for norm comparisons
      norm_tolerance = 1e-3
      same_norms = sumw_ref and all(abs(ao.sumW()/sumw_ref - 1) < norm_tolerance for ao in aos)
      ## Set the normto flag if possible (either because norms match or user-forced)
      if same_norms or opts.NORMALIZE_ALL:
          if not all(ao.annotations.has_key("ScaledBy") for ao in aos):
              print "WARNING: Abandoning normalized merge of path %s because not all inputs have ScaledBy attributes" % p
          else:
              ## Try to compute a target normalization from the 1/scalefactor-weighted norms of each run
              wtot = sum(1/float(ao.annotations["ScaledBy"]) for ao in aos)
              normto = sum(ao.sumW() / float(ao.annotations["ScaledBy"]) for ao in aos) / wtot
      ## Now that the normalization identifying heuristic is done, apply requested scaling
      scale = float(ao.annotation("yodamerge_scale"))
      if hasattr(ao, "scaleW") and scale != 1.0:
          ao.scaleW(scale)
      elif type(ao) is yoda.Scatter2D:
          ao.scale(1.0, scale)
      scale = ao.rmAnnotation("yodamerge_scale")

  ## Loop over and combine all data objects for this path
  output_object = None
  for ao in aos:
      ## Unscale first if normto != None
      # TODO: Move this into the normto calculation above?
      if normto:
          ao.scaleW(1/float(ao.annotations["ScaledBy"]))
      ## Combine as far as supported by this data type
      if output_object is None:
          output_object = ao
      elif hasattr(ao, "__iadd__"):
          output_object += ao
      else:
          print "WARNING: Analysis object %s of type %s cannot be merged" % (p, str(type(ao)))
          break
  ## Re-normalize after adding if normto != None
  if normto:
      output_object.normalize(normto)

  return output_object

parser = optparse.OptionParser(usage=__doc__)
parser.add_option('-o', '--output', default='-', dest='OUTPUT_FILE')
parser.add_option('-N', '--normalize-all', action="store_true", default=False, dest='NORMALIZE_ALL')
opts, fileargs = parser.parse_args()

## Put the incoming objects into a dict from each path to a list of histos and scalings
analysisobjects_in = {}
for fa in fileargs:
    filename, scale = fa, 1.0 
    if ":" in fa: 
        try:
            filename, scale = fa.rsplit(":", 1)
            scale = float(scale)
        except:
            sys.stderr.write("Error processing arg '%s' with file:scale format\n" % fa) 
    aos = yoda.read(filename)
    for aopath, ao in aos.iteritems():
        ao.setAnnotation("yodamerge_scale", scale)
        analysisobjects_in.setdefault(aopath, []).append(ao)

analysisobjects_out = {}
for path, aos in analysisobjects_in.iteritems():
  if(path in EFFICIENCIES): analysisobjects_out[path] = MergeEfficiency(path, aos)
  else: analysisobjects_out[path] = MergeDistribution(path, aos)

yoda.writeYODA(analysisobjects_out, opts.OUTPUT_FILE)
