# preView: A function for getting an exploratory preview of infection/exposure prevalence from sample data.
# Dependencies: dplyr, tidyr, ggplot2
# Inputs:
# --data..........a data.frame containing all the necessary data
# --outcome.......a character vector containing the (quoted) names of columns in `data` housing titers or exposure/infection status indicators (need not be binary)
# --by............a string containing the (quoted) name of the column in `data` containing the domain of the prevalence visualization (e.g., year or month); if a character vector of length >1, subsequent elements are appended to the head of `group`
# --is.interval...logical; indicates whether the domain (`by`) is an interval or categorical variable; determines whether lines and points or points only are used to plot prevalence
# --group.........a character vector containing the (quoted) names of columns in `data` housing additional groupings (e.g., sex, age); inherits extra categories from `by`; currently only 1 grouping is supported
# --alpha.........numeric; alpha level for calculating confidence intervals
# --method........character; determines the method used to calculate prevalence and confidence intervals (after Brown et al., 2001, Stat Sci 16:101--33).
# Outputs:
# --p.plot........returns a ggplot object that can be captured and edited

preView <- function( data, outcome, by, is.interval = TRUE, group = NULL, alpha = 0.05, method = c( "wald", "agresti-coull", "wilson", "jeffreys", "clopper-pearson" ) ){
  
  # Check to make sure the data has been supplied.
  if( missing( data ) ) warning( "Where's the data?" )
  
  # Check to see how many domains (`by` groupings) were supplied; pass extras on to `group` if necessary.
  if( ( n.by <- length( by ) ) > 1 ){
    warning( "Only one domain is supported; transferring secondaries to group..." )
    group <- c( by[ 2:n.by ], group )
    by <- by[ 1 ]
    n.by <- 1
  }
  
  # Check to see how many `group` groupings were supplied (or were inherited from `by`); ignore extras if they exist.
  if( ( n.groups <- length( group ) ) > 1 ){
    warning( "Only one grouping supported; picking the first given..." )
    group <- group[ 1 ]
    n.groups <- 1
  }
  
  ##### MAKE CALCULATIONS #####
  
  # calculate kappa, the z-score defined as the inverse of the Gaussian distribution function (quantile function) at 1-(a/2) (For a=0.05, k=1.959...)
  k <- qnorm( p = ( 1 - ( alpha / 2 ) ) )
  
  # Identify the method used to calculate confidence intervals (and centers)--following Brown et al. (2001)
  method <- match.arg( method )
  
  # Define functions used to calculate prevalence/proportions by different `method`s.
  p.fun <- switch( method, 
                   "wald" = function( .t ){ sum( as.logical( .t ), na.rm = TRUE ) / sum( !is.na( .t ) ) }, 
                   "agresti-coull" = function( .t ){ ( sum( as.logical( .t ), na.rm = TRUE ) + ( 0.5 * k^2 ) ) / ( sum( !is.na( .t ) ) + k^2 ) }, 
                   "wilson" = function( .t ){ ( sum( as.logical( .t ), na.rm = TRUE ) + ( 0.5 * k^2 ) ) / ( sum( !is.na( .t ) ) + k^2 ) }, 
                   "jeffreys" = function( .t ){ sum( as.logical( .t ), na.rm = TRUE ) / sum( !is.na( .t ) ) }, 
                   "clopper-pearson" = function( .t ){ sum( as.logical( .t ), na.rm = TRUE ) / sum( !is.na( .t ) ) } )
  
  # Define functions used to calculate confidence intervals by different `method`s.
  ci.fun <- switch( method, 
                    "wald" = function( .x, .p, .k, .n, .l = c( "l", "u" ) ){ 
                      if( identical( .l, "l" ) )
                        .p - ( .k * sqrt( ( .p * ( 1 - .p ) ) / .n ) )
                      else .p + ( .k * sqrt( ( .p * ( 1 - .p ) ) / .n ) ) }, 
                    "agresti-coull" = function( .p, .k, .n, .l = c( "l", "u" ) ){ 
                      if( identical( .l, "l" ) )
                        .p - ( .k * sqrt( ( .p * ( 1 - .p ) ) / .n ) ) 
                      else .p + ( .k * sqrt( ( .p * ( 1 - .p ) ) / .n ) ) },
                    "wilson" = function( .x, .p, .k, .n, .l = c( "l", "u" ) ){ 
                      if( identical( .l, "l" ) ) 
                        .p - ( ( ( .k * sqrt( .n ) ) / ( .n + .k^2 ) ) * sqrt( ( .p * ( 1 - .p ) ) + ( .k^2 / ( 4 * .n ) ) ) )
                      else .p + ( ( ( .k * sqrt( .n ) ) / ( .n + .k^2 ) ) * sqrt( ( .p * ( 1 - .p ) ) + ( .k^2 / ( 4 * .n ) ) ) ) }, 
                    "jeffreys" = function( .x, .p, .k, .n, .l = c( "l", "u" ) ){
                      if( identical( .l, "l" ) )
                        qbeta( alpha / 2, .x + 0.5, .n - .x + 0.5 )
                      else qbeta( 1 - ( alpha / 2 ), .x + 0.5, .n - .x + 0.5 ) }, 
                    "clopper-pearson" = function( .x, .p, .k, .n, .l = c( "l", "u" ) ){ 
                      if( identical( .l, "l" ) )
                        ifelse( .x == 0, 0, qbeta( alpha / 2, .x, .n - .x + 1 ) )
                      else ifelse( .x == .n, 1, qbeta( 1 - ( alpha / 2 ), .x + 1, .n - .x ) ) } )
  
  # Summarize the data.
  p.summary <- data %>%                                                             # From the data,
    dplyr::select_( .dots = c( outcome, by, group ) ) %>%                           # Filter out extraneous variables.
    dplyr::ungroup() %>%                                                            # Make sure there aren't any residual groupings.
    dplyr::group_by_( .dots = c( by, group ) ) %>%                                  # Then establish groupings by `by` and `group`.
    dplyr::summarize_at( .cols = outcome,                                           # For each outcome,
                         .funs = funs( x = sum( as.logical( . ), na.rm = TRUE ),    # calculate the number infected or exposed (as.logical converts to binary),
                                       n = sum( !is.na( . ) ),                      # the total number sampled,
                                       p = p.fun( . ) ) )                           # and the prevalence (varies by `method`; for "wald", a simple proportion of x/n)
  
  if( ( n.outcome <- length( outcome ) ) == 1 ){                                         # dplyr::summarize creates column headings differently
    p.summary <- p.summary %>%                                                           # depending on whether or not it is performed at multiple .cols
      dplyr::rename_( .dots = c( setNames( "x", paste( outcome, "x", sep = "_" ) ),      # if only one `outcome` is specified, I need to adjust the column headings
                                 setNames( "n", paste( outcome, "n", sep = "_" ) ),      # to prepare the summary table for the next step.
                                 setNames( "p", paste( outcome, "p", sep = "_" ) ) ) )   # Essentially, statistic becomes outcome_statistic.
  }
  
  p.summary <- p.summary %>%                                                             # Transform the preliminary data summary
    tidyr::gather( key = disease_parameter,                                              # into long-format, gathering the outcome_statistic columns
                   value = value, 
                   -dplyr::one_of( c( by, group ) ) ) %>% 
    tidyr::separate( col = disease_parameter,                                            # Now separate `outcome` from `statistic, making two columns,
                     into = c( "disease", "parameter" ), 
                     sep = "_" ) %>% 
    tidyr::spread( key = parameter,                                                      # and spread it back out, keeping `outcome` in a column for plot grouping.
                   value = value ) %>% 
    dplyr::mutate( ci.hi = ci.fun( x, p, k, n, "u" ),                                      # Calculate CIs based on the appropriate `method`ology.
                   ci.lo = ci.fun( x, p, k, n, "l" ) )
  
  ##### MAKE PLOT #####
  
  if( !is.null( group ) ){                                                               # First define the aesthetic mapping:
    mapping <- ggplot2::aes_( x = as.name( by ), 
                              y = quote( p ), 
                              shape = as.name( group ),                                  # If a `group` is specified, it gets mapped to shape (and linetype, if is.interval=TRUE).
                              colour = quote( disease ) )
  } else {
    mapping <- ggplot2::aes_( x = as.name( by ),                                         # Either way, `by` gets mapped to the x-domain,
                              y = quote( p ),                                            # prevalence is plotted on the y-axis,
                              colour = quote( disease ) )                                # and `outcome`s are distinguished by color.
  }
  
  p.plot <- ggplot2::ggplot( data = p.summary, 
                             mapping = mapping ) + 
    ggplot2::geom_linerange( linetype = "solid",                                         # A linerange is used to indicate confidence intervals,
                             ggplot2::aes( ymin = ci.lo, 
                                           ymax = ci.hi ), 
                             position = ggplot2::position_dodge( width = 0.3 ) ) +       # with position dodging to prevent overplotting.
    ggplot2::geom_point( position = ggplot2::position_dodge( width = 0.3 ) ) +           # Points indicate mean prevalence/proportion (again, with dodging).
    ggplot2::theme_classic()
  
  if( is.interval ){                                                                     # If `by` is interval, lines connect sequential observations, with linetype mapped to `group` (if specified).
    p.plot <- p.plot + ggplot2::geom_line( ggplot2::aes_( linetype = if( !is.null( group ) ) as.name( group ) else NULL ), 
                                           position = ggplot2::position_dodge( width = 0.3 ) )
  }
  
  list( p.data = p.summary,                                                              # Return a list containing the summarized data as well as the plot.
        p.plot = p.plot )
  
}