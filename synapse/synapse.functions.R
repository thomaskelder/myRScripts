############################################
## Utility functions for Synapse R client ##
############################################
createOrGetEntity = function(entity) {
  name = properties(entity)$name
  id = synapseQuery(paste("select id from entity where name == '", name, "'", sep=''))
  
  if(is.null(id)) { #Doesn't exist yet, create
    entity = createEntity(entity)
  } else { #Already exists, get
    entity = getEntity(as.character(id))
  }
  entity
}

filesFromEntity = function(entity) {
  paste(entity$cacheDir, entity$files, sep="/")
}

loadSynapseData = function(id) {
  ent = downloadEntity(id)
  load(filesFromEntity(ent), .GlobalEnv)
}