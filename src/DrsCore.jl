module DrsCore

using Logging
@Logging.configure(level=DEBUG)

export CHUZR

function CHUZR()
	@debug("CHUZR")
end

end
